/*
 * Binning with MetaBAT2, MaxBin2 and CONCOCT
 */
include { INDEX_FASTA                 } from '../../modules/local/bwa_mem2/index'
include { FEATURED_ALIGNMENT          } from '../../modules/local/bwa_mem2/binning_alignment'
include { METABAT2_METABAT2           } from '../../modules/nf-core/metabat2/metabat2/main'
include { MAXBIN2                     } from '../../modules/nf-core/maxbin2/main'
include { CONVERT_DEPTHS              } from '../../modules/local/mag/convert_depths'
include { CONCOCT_SUBWF               } from './concoct-subwf'

/*
 * Subworkflow
 */

workflow BINNING {

    take:
    run_concoct  // channel: [ val(meta), path(assembly), [reads]]
    run_metabat
    run_maxbin
    run_depth

    main:

    ch_versions = Channel.empty()

    // filter nulls
    input_concoct_bins = run_concoct.filter { it != null }
    input_maxbin_bins = run_maxbin.filter { it != null }
    input_metabat_bins = run_metabat.filter { it != null }
    input_depth = run_depth.filter { it != null }

    // Create a unified channel of all samples that need binning
    // This ensures each unique sample is only processed once for common steps
    all_binning_samples = input_concoct_bins
        .mix(input_maxbin_bins)
        .mix(input_metabat_bins)
        .mix(input_depth)
        .unique { meta, _assembly, _reads -> meta.id } // Remove duplicates based on sample ID

    /*
    // Run common pre-binning steps for all samples that require any binning
    */
    INDEX_FASTA(
        all_binning_samples.map { meta, assembly, _reads -> [meta, assembly] }
    )
    ch_versions = ch_versions.mix( INDEX_FASTA.out.versions )

    /*
    * ---
    * optional args: run_jgi_depth = true
    *                run_concoct_table_generation = true
    *                concoct_generate_bed_file = true
    * output: depth ( produced by jgi_summarize_bam_contig_depths )
    *         TSV from concoct_coverage_table.py, FASTA from cut_up_fasta.py
    *         idxstats for multiqc
    */
    FEATURED_ALIGNMENT(
        all_binning_samples.map { meta, _assembly, reads -> [meta, reads] }.join( INDEX_FASTA.out.fasta_with_index ),
        true,
        true,
        true
    )
    ch_versions = ch_versions.mix( FEATURED_ALIGNMENT.out.versions )

    /*
    * --- binning with CONCOCT ---
    * CONCOCT -> merge clusters -> extract fasta bins
    * input: [ val(meta), concoct_tsv, concoct_fasta, assembly_fasta ]
    */
    CONCOCT_SUBWF(
       FEATURED_ALIGNMENT.out.concoct_data
        .join(all_binning_samples.map { meta, assembly, _reads -> [meta, assembly] })
        .join(input_concoct_bins)
        .map{ meta, tsv, concoct_fasta, uncompressed_assembly, _assembly, _reads -> [meta, tsv, concoct_fasta, uncompressed_assembly]}
    )
    ch_versions = ch_versions.mix( CONCOCT_SUBWF.out.versions )

    concoct_output = CONCOCT_SUBWF.out.bins

    /*
    * --- binning with METABAT2 ---
    */
    METABAT2_METABAT2 (
        input_metabat_bins.map { meta, assembly, _reads -> [meta, assembly] }.join(FEATURED_ALIGNMENT.out.depth)
    )
    ch_versions = ch_versions.mix( METABAT2_METABAT2.out.versions )

    metabat_output = METABAT2_METABAT2.out.fasta

    /* 
    * --- binning with MAXBIN2 ---
    * only used for prokaryotic MAGs
    */
    maxbin_output  = Channel.empty()

    if ( !params.skip_prok ) {
        // convert metabat2 depth files to maxbin2
        CONVERT_DEPTHS (
            input_maxbin_bins.map { meta, assembly, _reads -> [meta, assembly] }.join(FEATURED_ALIGNMENT.out.depth)
        )
        ch_versions = ch_versions.mix( CONVERT_DEPTHS.out.versions )

        /*
        * maxbin process
        * we provide empty reads as we don't want maxbin2 to calculate for us.
        */
        MAXBIN2 (
           CONVERT_DEPTHS.out.output.map { meta, assembly, depth -> [meta, assembly, [], depth ] }
        )  // output can be empty folder
        ch_versions = ch_versions.mix( MAXBIN2.out.versions )

        maxbin_output = MAXBIN2.out.binned_fastas
    }
    
    emit:
    maxbin_bins      = maxbin_output
    concoct_bins     = concoct_output
    metabat_bins     = metabat_output
    jgi_depth        = FEATURED_ALIGNMENT.out.depth
    versions         = ch_versions
}