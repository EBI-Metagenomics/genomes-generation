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

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { FINALIZE_LOGGING            } from '../modules/local/utils'
include { GUNZIP as GUNZIP_ASSEMBLY   } from '../modules/local/utils'
//include { MULTIQC                     } from '../modules/nf-core/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BINNING              } from '../subworkflows/local/binning'
include { DECONTAMINATION      } from '../subworkflows/local/decontamination'
//include { EUKCC_MERGE          } from '../subworkflows/local/eukcc_merge'
include { EUK_MAGS_GENERATION  } from './eukaryotic_mags_generation'
include { INPUT_PREPROCESSING  } from '../subworkflows/local/input_preprocessing'
//include { PREPARE_UPLOAD_FILES } from '../subworkflows/local/prepare_upload'
//include { PROK_MAGS_GENERATION } from '../subworkflows/local/prok_mags_generation'
include { QC_AND_MERGE_READS   } from '../subworkflows/local/qc_and_merge'

/*
~~~~~~~~~~~~~~~~~~
    DBs
~~~~~~~~~~~~~~~~~~
*/
//cat_db_folder      = file(params.cat_db_folder, checkIfExists: true)
//cat_diamond_db     = file(params.cat_diamond_db, checkIfExists: true)
//cat_taxonomy_db    = file(params.cat_taxonomy_db, checkIfExists: true)
//gunc_db            = file(params.gunc_db, checkIfExists: true)
//checkm2_db         = file(params.checkm2_db, checkIfExists: true)
//busco_db           = file(params.busco_db, checkIfExists: true)
//gtdbtk_db          = file(params.gtdbtk_db, checkIfExists: true)
//rfam_rrna_models   = file(params.rfam_rrna_models, checkIfExists: true)

//assembly_software  = file(params.assembly_software_file, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Run workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GGP {

    take:
    samplesheet

    main:

    ch_versions    = Channel.empty()
    //ch_log         = Channel.empty()

    //euk_genomes    = Channel.empty()
    //stats_euks     = Channel.empty()
    //coverage_euks  = Channel.empty()
    //prok_genomes   = Channel.empty()
    //stats_proks    = Channel.empty()
    //coverage_proks = Channel.empty()
    //rna            = Channel.empty()
    //taxonomy_euks  = Channel.empty()
    //taxonomy_proks = Channel.empty()

    /*
    * ---- combine data for reads, contigs and bins ----
    */
    samplesheet.multiMap{ meta, fq1, fq2, assembly, concoct, metabat, maxbin, depth ->
        def is_single_end = (fq2 == [])
        def reads = is_single_end ? [fq1] : [fq1, fq2]
        assembly_and_reads : tuple(meta + [single_end: is_single_end], assembly, reads)
        concoct: concoct ? tuple(meta + [single_end: is_single_end], concoct) : null
        metabat: metabat ? tuple(meta + [single_end: is_single_end], metabat) : null
        maxbin: maxbin ? tuple(meta + [single_end: is_single_end], maxbin) : null
        jgi_depth: depth ? tuple(meta + [single_end: is_single_end], depth) : null
    }.set {
        input
    }
    concoct_sample_ids = input.concoct.filter { it != null }
        .map { meta, f -> meta.id }
        .collect()
        .ifEmpty([])
    metabat_sample_ids = input.metabat.filter { it != null }
        .map { meta, f -> meta.id }
        .collect()
        .ifEmpty([])
    maxbin_sample_ids = input.maxbin.filter { it != null }
        .map { meta, f -> meta.id }
        .collect()
        .ifEmpty([])
    depth_sample_ids = input.jgi_depth.filter { it != null }
        .map { meta, d -> meta.id }
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
        input.assembly_and_reads
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

    // --- make data strucrute for binning
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
    //branched_for_tools.run_maxbin.filter { it != null }

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

    existing_concoct = branched_for_tools.include_concoct.join(input.concoct).filter { it != null }.map{ meta, _1, _2, bins -> [meta, bins] }
    existing_maxbin = branched_for_tools.include_maxbin.join(input.maxbin).filter { it != null }.map{ meta, _1, _2, bins -> [meta, bins] }
    existing_metabat = branched_for_tools.include_metabat.join(input.metabat).filter { it != null }.map{ meta, _1, _2, bins -> [meta, bins] }
    existing_depth = branched_for_tools.include_depth.join(input.jgi_depth).filter { it != null }.map{ meta, _1, _2, depth -> [meta, depth] }

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
        ch_versions = ch_versions.mix( EUK_MAGS_GENERATION.out.versions )
    }
    // for multiqc
    // binning samtools for multiqc
    // TODO return fastqc from INPUT_PREPROCESSING
    // TODO return fastqc result after decontamination
    // TODO return ALIGN samtools stats

}
