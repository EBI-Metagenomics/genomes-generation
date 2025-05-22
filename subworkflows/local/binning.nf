/*
 * Binning with MetaBAT2, MaxBin2 and CONCOCT
 */
include { INDEX_FASTA                 } from '../../modules/local/bwa_mem2/index'
include { FEATURED_ALIGNMENT          } from '../../modules/local/bwa_mem2/binning_alignment'
include { GUNZIP as GUNZIP_ASSEMBLY   } from '../../modules/local/utils'
include { METABAT2_METABAT2           } from '../../modules/nf-core/metabat2/metabat2/main'
include { MAXBIN2                     } from '../../modules/nf-core/maxbin2/main'
include { CONVERT_DEPTHS              } from '../../modules/local/mag/convert_depths'
include { CONCOCT_SUBWF               } from './concoct-subwf'

/*
 * Subworkflow
 */

workflow BINNING {

    take:
    assembly_and_reads  // channel: [ val(meta), path(assembly), [reads] ]

    main:

    ch_versions = Channel.empty()
    // optional //
    metabat_output = Channel.empty()
    concoct_output = Channel.empty()
    maxbin_output  = Channel.empty()

    /*
    * --- create index for assembly.fasta ---
    */
    INDEX_FASTA( assembly_and_reads.map { meta, assembly, reads -> [meta, assembly]} )
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
        assembly_and_reads.map { meta, assembly, reads -> [ meta, reads ] }.join( INDEX_FASTA.out.fasta_with_index ),
        true,
        true,
        true
    )
    ch_versions = ch_versions.mix( FEATURED_ALIGNMENT.out.versions )

    /*
    * --- uncompress assembly fasta ---
    */
    GUNZIP_ASSEMBLY(
        assembly_and_reads.map { meta, assembly, _ -> [ meta, assembly ] }
    )

    /*
    * --- binning with MAXBIN2 ---
    */
    if ( !params.skip_maxbin2 ) {
        // convert metabat2 depth files to maxbin2
        CONVERT_DEPTHS (
            assembly_and_reads.map { meta, assembly, _ -> [ meta, assembly ] }.join(FEATURED_ALIGNMENT.out.depth)
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

    /*
    * --- binning with METABAT2 ---
    */
    if ( !params.skip_metabat2 ) {

        METABAT2_METABAT2 (
            assembly_and_reads.map { meta, assembly, _ -> [ meta, assembly ] }.join(FEATURED_ALIGNMENT.out.depth)
        )
        ch_versions = ch_versions.mix( METABAT2_METABAT2.out.versions )

        metabat_output = METABAT2_METABAT2.out.fasta
    }

    /*
    * --- binning with CONCOCT ---
    */
    if ( !params.skip_concoct ) {
        /*
        * CONCOCT -> merge clusters -> extract fasta bins
        * input: [ val(meta), concoct_tsv, concoct_fasta, assembly_fasta ]
        */
        CONCOCT_SUBWF(
           FEATURED_ALIGNMENT.out.concoct_data.join(assembly_and_reads.map { meta, assembly, _ -> [ meta, assembly ] })
        )
        ch_versions = ch_versions.mix( CONCOCT_SUBWF.out.versions )

        concoct_output = CONCOCT_SUBWF.out.bins
    }

    emit:
    maxbin_bins      = maxbin_output
    concoct_bins     = concoct_output
    metabat_bins     = metabat_output
    versions         = ch_versions
}