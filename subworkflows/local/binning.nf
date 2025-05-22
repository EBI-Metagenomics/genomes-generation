/*
 * Binning with MetaBAT2, MaxBin2 and CONCOCT
 */
include { INDEX_FASTA        } from '../../modules/local/bwa_mem2/index'
include { FEATURED_ALIGNMENT } from '../../modules/local/bwa_mem2/binning_alignment'
include { METABAT2_METABAT2  } from '../../modules/nf-core/metabat2/metabat2/main'
include { MAXBIN2            } from '../../modules/nf-core/maxbin2/main'
include { CONVERT_DEPTHS     } from '../../modules/local/mag/convert_depths'
include { CONCOCT_SUBWF      } from './concoct-subwf'

workflow BINNING {

    take:
    assembly_and_reads  // channel: [ val(meta), path(assembly), depth, concoct_tsv, concoct_fasta ]

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
    */
    FEATURED_ALIGNMENT(
        assembly_and_reads.map { meta, assembly, reads -> [ meta, reads ] }.join( INDEX_FASTA.out.fasta_with_index ),
        true,
        true,
        true
    )
    ch_versions = ch_versions.mix( FEATURED_ALIGNMENT.out.versions )


    assemblies_depth_coverage.multiMap { meta, assembly, depth, coverage, concoct_fasta ->
        assembly: [ meta, assembly ]
        concoct_tsv: [ meta, coverage, concoct_fasta ]
        depth: [meta, depth]
    }.set {
        input
    }

    // convert metabat2 depth files to maxbin2
    if ( !params.skip_maxbin2 ) {

        CONVERT_DEPTHS ( input.assembly.join(input.depth) )

        ch_maxbin2_input = CONVERT_DEPTHS.out.output.map { meta, assembly, depth ->
            // we provide empty reads as we don't want maxbin2 to calculate for us.
            [meta, assembly, [], depth ]
        }

        MAXBIN2 ( ch_maxbin2_input )  // output can be empty folder
        maxbin_output = MAXBIN2.out.binned_fastas

        ch_versions = ch_versions.mix( CONVERT_DEPTHS.out.versions.first() )
        ch_versions = ch_versions.mix( MAXBIN2.out.versions.first() )
    }

    if ( !params.skip_metabat2 ) {

        METABAT2_METABAT2 ( input.assembly.join(input.depth) )

        metabat_output = METABAT2_METABAT2.out.fasta

        ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())
    }

    if ( !params.skip_concoct ) {

        CONCOCT_SUBWF( input.concoct_tsv.join(input.assembly) )

        concoct_output = CONCOCT_SUBWF.out.bins

        ch_versions = ch_versions.mix( CONCOCT_SUBWF.out.versions.first() )
    }

    emit:
    maxbin_bins      = maxbin_output
    concoct_bins     = concoct_output
    metabat_bins     = metabat_output
    versions         = ch_versions
}