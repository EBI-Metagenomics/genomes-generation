/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap    } from 'plugin/nf-schema'  // for multiqc

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { FINALIZE_LOGGING            } from '../modules/local/utils'
include { PIGZ as COMPRESS_BINS       } from '../modules/local/compress/pigz'
//include { MULTIQC                     } from '../modules/nf-core/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EUK_MAGS_GENERATION        } from './eukaryotic_mags_generation'
include { PROK_MAGS_GENERATION       } from './prokaryotic_mags_generation'

include { BINNING                    } from '../subworkflows/local/binning'
include { DECONTAMINATION            } from '../subworkflows/local/decontamination'
include { INPUT_PREPROCESSING        } from '../subworkflows/local/input_preprocessing'
include { QC_AND_MERGE_READS         } from '../subworkflows/local/qc_and_merge'
include { UPLOAD_MAGS                } from '../subworkflows/local/upload'
include { UPLOAD_MAGS as UPLOAD_BINS } from '../subworkflows/local/upload'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Run workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GGP {

    take:
    assembly_and_reads
    input_binners_concoct
    input_binners_metabat
    input_binners_maxbin
    input_jgi_depth
    assembly_software_file

    main:

    ch_versions    = Channel.empty()
    ch_log         = Channel.empty()

    euk_mags       = Channel.empty()
    euk_bins       = Channel.empty()
    stats_euks     = Channel.empty()
    coverage_euks  = Channel.empty()
    taxonomy_euks  = Channel.empty()
    prok_mags      = Channel.empty()
    prok_bins      = Channel.empty()
    stats_proks    = Channel.empty()
    coverage_proks = Channel.empty()
    rna_proks      = Channel.empty()
    taxonomy_proks = Channel.empty()

    /*
    * ---- Optional input binners parsing ----
    */
    concoct_sample_ids = input_binners_concoct.filter { it != null }
        .map { meta, _files -> meta.id }
        .collect()
        .ifEmpty([])
    metabat_sample_ids = input_binners_metabat.filter { it != null }
        .map { meta, _files -> meta.id }
        .collect()
        .ifEmpty([])
    maxbin_sample_ids = input_binners_maxbin.filter { it != null }
        .map { meta, _files -> meta.id }
        .collect()
        .ifEmpty([])
    depth_sample_ids = input_jgi_depth.filter { it != null }
        .map { meta, _depth -> meta.id }
        .collect()
        .ifEmpty([])

    // -------------------- WORKFLOW --------------------
    /*
    * --- pre-processing input files ---
    * skip that step with --skip_preprocessing_input
    * change ERR to ERZ in reads
    * change . to _ in assembly files
    */
    INPUT_PREPROCESSING(
        assembly_and_reads
    )
    ch_versions = ch_versions.mix(INPUT_PREPROCESSING.out.versions)

    /*
    * --- trimming reads ----
    * merge step is regulated with --merge_pairs (default: false)
    */
    QC_AND_MERGE_READS(
        INPUT_PREPROCESSING.out.assembly_and_reads.map{ meta, _1, reads -> [meta, reads] }
    )
    ch_versions = ch_versions.mix( QC_AND_MERGE_READS.out.versions )

    /*
    * --- reads decontamination ----
    * skip that step with --skip_decontamination
    */
    reference_genome         = file(params.ref_genome, checkIfExists: true)
    reference_genome_index   = file("${reference_genome.parent}/*.fa*.*")
    DECONTAMINATION(
        QC_AND_MERGE_READS.out.reads,
        reference_genome,
        reference_genome_index
    )
    ch_versions = ch_versions.mix( DECONTAMINATION.out.versions )

    // --- check filtered reads quality --- //
    FASTQC (
        DECONTAMINATION.out.decontaminated_reads
    )
    ch_versions = ch_versions.mix( FASTQC.out.versions )

    // --- make data structure for binning
    tool_availability = INPUT_PREPROCESSING.out.assembly_and_reads.map{ meta, assembly, _2 -> [meta, assembly] }
        .join(DECONTAMINATION.out.decontaminated_reads)
        .combine(concoct_sample_ids.map { [it] })
        .combine(metabat_sample_ids.map { [it] })
        .combine(maxbin_sample_ids.map { [it] })
        .combine(depth_sample_ids.map { [it] })
        .map { meta, assembly, reads, concoct_list, metabat_list, maxbin_list, depth_list ->
            def tools_presented = []
            if (concoct_list.contains(meta.id)) tools_presented.add('concoct')
            if (metabat_list.contains(meta.id)) tools_presented.add('metabat')
            if (maxbin_list.contains(meta.id)) tools_presented.add('maxbin')
            if (depth_list.contains(meta.id)) tools_presented.add('depth')
            return [meta, assembly, reads, tools_presented]
        }

    tool_availability
        .multiMap { meta, assembly, reads, tools ->
            run_concoct: !tools.contains('concoct') ? [meta, assembly, reads] : null
            include_concoct: tools.contains('concoct') ? [meta, assembly, reads] : null
            run_metabat: !tools.contains('metabat') ? [meta, assembly, reads] : null  
            include_metabat: tools.contains('metabat') ? [meta, assembly, reads] : null  
            run_maxbin: !tools.contains('maxbin') ? [meta, assembly, reads] : null
            include_maxbin: tools.contains('maxbin') ? [meta, assembly, reads] : null
            run_depth: !tools.contains('depth') ? [meta, assembly, reads] : null
            include_depth: tools.contains('depth') ? [meta, assembly, reads] : null
        }
        .set { branched_for_tools }

    /*
    * --- binning ---
    */
    BINNING(
        branched_for_tools.run_concoct,
        branched_for_tools.run_metabat,
        branched_for_tools.run_maxbin,
        branched_for_tools.run_depth,
    )
    ch_versions = ch_versions.mix( BINNING.out.versions )

    existing_concoct = branched_for_tools.include_concoct.join(input_binners_concoct).filter { it != null }.map{ meta, _1, _2, bins -> [meta, bins] }
    existing_maxbin = branched_for_tools.include_maxbin.join(input_binners_maxbin).filter { it != null }.map{ meta, _1, _2, bins -> [meta, bins] }
    existing_metabat = branched_for_tools.include_metabat.join(input_binners_metabat).filter { it != null }.map{ meta, _1, _2, bins -> [meta, bins] }
    existing_depth = branched_for_tools.include_depth.join(input_jgi_depth).filter { it != null }.map{ meta, _1, _2, depth -> [meta, depth] }

    all_concoct_bins = existing_concoct.mix(BINNING.out.concoct_bins)
    all_maxbin_bins = existing_maxbin.mix(BINNING.out.maxbin_bins)
    all_metabat_bins = existing_metabat.mix(BINNING.out.metabat_bins)
    all_jgi_depth = existing_depth.mix(BINNING.out.jgi_depth).groupTuple().map{ meta, depth_list -> [ meta, depth_list[0] ] }

    /*
    * --- MAGs generation ---
    */ 
    if ( !params.skip_euk ) {
        EUK_MAGS_GENERATION(
            INPUT_PREPROCESSING.out.assembly_and_reads.map{ meta, assembly, _2 -> [meta, assembly] }
            .join(DECONTAMINATION.out.decontaminated_reads)
            .join(all_concoct_bins)
            .join(all_metabat_bins)
            .join(all_jgi_depth)
        )
        euk_mags      = euk_mags.mix( EUK_MAGS_GENERATION.out.mags_fastas )
        euk_bins      = euk_bins.mix( EUK_MAGS_GENERATION.out.bins_fastas )
        stats_euks    = stats_euks.mix( EUK_MAGS_GENERATION.out.stats )
        coverage_euks = coverage_euks.mix( EUK_MAGS_GENERATION.out.coverage )
        taxonomy_euks = taxonomy_euks.mix( EUK_MAGS_GENERATION.out.taxonomy )
        ch_versions   = ch_versions.mix( EUK_MAGS_GENERATION.out.versions )
        ch_log        = ch_log.mix( EUK_MAGS_GENERATION.out.progress_log )
    }

    if ( !params.skip_prok ) {
        // input: tuple( meta, concoct, metabat, maxbin, depth_file)
        collected_binners_assembly_and_depth = all_concoct_bins
            .join( all_maxbin_bins, remainder: true ) 
            .join( all_metabat_bins, remainder: true ) 
            .join( INPUT_PREPROCESSING.out.assembly_and_reads.map{ meta, assembly, _reads -> [meta, assembly] }, remainder: true ) 
            .join( all_jgi_depth, remainder: true )

        PROK_MAGS_GENERATION(
            collected_binners_assembly_and_depth
        )

        prok_mags      = prok_mags.mix( PROK_MAGS_GENERATION.out.mags_fastas )
        prok_bins      = prok_bins.mix( PROK_MAGS_GENERATION.out.bins_fastas )
        stats_proks    = stats_proks.mix( PROK_MAGS_GENERATION.out.stats )
        coverage_proks = coverage_proks.mix( PROK_MAGS_GENERATION.out.coverage )
        rna_proks      = rna_proks.mix( PROK_MAGS_GENERATION.out.mags_rna )
        taxonomy_proks = taxonomy_proks.mix( PROK_MAGS_GENERATION.out.taxonomy )
        ch_versions    = ch_versions.mix( PROK_MAGS_GENERATION.out.versions )
        ch_log         = ch_log.mix( PROK_MAGS_GENERATION.out.progress_log )
    }
    
    if ( params.skip_euk && params.skip_prok ) {
        log.error "You skipped both prokaryotes and eukaryotes. No results for MAGs. Exiting."
        exit(1)
    }

    if ( params.upload_mags ) {
        UPLOAD_MAGS (
            euk_mags.ifEmpty([]),
            prok_mags.ifEmpty([]),
            assembly_software_file,
            stats_euks.ifEmpty([]),
            stats_proks.ifEmpty([]),
            coverage_euks.ifEmpty([]),
            coverage_proks.ifEmpty([]),
            rna_proks.ifEmpty([]),
            taxonomy_euks.ifEmpty([]),
            taxonomy_proks.ifEmpty([]),
            "mags")  // mags_or_bins_flag = "mags"
        ch_versions = ch_versions.mix( UPLOAD_MAGS.out.versions )
    }

    if ( params.upload_bins ) {
        UPLOAD_BINS (
                euk_bins.ifEmpty([]),
                prok_bins.ifEmpty([]),
                assembly_software_file,
                stats_euks.ifEmpty([]),
                stats_proks.ifEmpty([]),
                coverage_euks.ifEmpty([]),
                coverage_proks.ifEmpty([]),
                [],   // empty input, because we skip rRNA prediction for bins
                taxonomy_euks.ifEmpty([]),
                taxonomy_proks.ifEmpty([])
                , "bins")  // mags_or_bins_flag = "bins"
        ch_versions = ch_versions.mix( UPLOAD_BINS.out.versions )
    }


    // for multiqc
    // binning samtools for multiqc
    // TODO return fastqc from INPUT_PREPROCESSING
    // TODO return fastqc result after decontamination
    // TODO return ALIGN samtools stats

    emit:
    versions          = ch_versions        // channel: [ versions.yml ]
    pipeline_logging  = ch_log.collectFile(name: 'pipeline_logging.txt')
}
