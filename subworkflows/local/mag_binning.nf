/*
 * Binning with MetaBAT2, MaxBin2 and Concoct
 */
include { METABAT2_METABAT2                     } from '../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS  } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { MAXBIN2                               } from '../../modules/nf-core/maxbin2/main'

include { CONVERT_DEPTHS                        } from '../../modules/local/mag/convert_depths'
include { RENAME_MAXBIN                         } from '../../modules/local/rename_maxbin/main'
include { FASTA_BINNING_CONCOCT                 } from '../nf-core/fasta_binning_concoct/main'

workflow BINNING {

    take:
    assemblies_bams  // channel: [ val(meta), path(assembly), path(bams), path(bais) ]

    main:

    ch_versions = Channel.empty()
    // optional //
    metabat_output = Channel.empty()
    concoct_output = Channel.empty()

    // generate coverage depths for each contig
    ch_summarizedepth_input = assemblies_bams.map { meta, assembly, bams, bais ->
        [ meta, bams, bais ]
    }

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )

    ch_metabat_depths = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            [ meta + [binner: 'MetaBAT2'], depths ]
        }

    ch_versions = ch_versions.mix( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first() )

    // combine depths back with assemblies
    ch_metabat2_input = assemblies_bams
        .map { meta, assembly, bams, bais ->
            [meta + [binner: 'MetaBAT2'], assembly, bams, bais]
        }
        .join( ch_metabat_depths )
        .map { meta, assembly, bams, bais, depths ->
            [ meta, assembly, depths ]
        }

    // convert metabat2 depth files to maxbin2
    if ( !params.skip_maxbin2 ) {
        
        CONVERT_DEPTHS ( ch_metabat2_input )

        ch_maxbin2_input = CONVERT_DEPTHS.out.output.map { meta, assembly, depth ->
            // we provide empty reads as we don't want maxbin2 to calculate for us.
            [meta.subMap('id', 'erz') + [binner: 'MaxBin2'], assembly, [], depth ]
        }

        MAXBIN2 ( ch_maxbin2_input )

        RENAME_MAXBIN ( MAXBIN2.out.binned_fastas )

        maxbin_output = RENAME_MAXBIN.out.renamed_bins.map { meta, bins ->
            [meta.subMap('id', 'erz'), bins]
        }

        ch_versions = ch_versions.mix( CONVERT_DEPTHS.out.versions.first() )
        ch_versions = ch_versions.mix( MAXBIN2.out.versions.first() )
        ch_versions = ch_versions.mix( RENAME_MAXBIN.out.versions.first() )
    }

    if ( !params.skip_metabat2 ) {

        METABAT2_METABAT2 ( ch_metabat2_input )

        metabat_output = METABAT2_METABAT2.out.fasta.map { meta, bins ->
            [ meta.subMap('id', 'erz'), bins ]
        }

        ch_versions.mix(METABAT2_METABAT2.out.versions.first())
    }

    if ( !params.skip_concoct ){

        assemblies_bams.map { meta, assembly, bams, bais ->
            [ meta + [binner: 'CONCOCT'], assembly, bams, bais ]
        }.multiMap { meta, assembly, bams, bais ->
            bins: [ meta, assembly ]
            bams: [ meta, bams, bais ]
        }.set { ch_concoct_input }

        FASTA_BINNING_CONCOCT( ch_concoct_input.bins, ch_concoct_input.bams )

        concoct_output = FASTA_BINNING_CONCOCT.out.bins.map { meta, bins ->
            [ meta.subMap('id', 'erz'), bins ]
        }

        ch_versions = ch_versions.mix( FASTA_BINNING_CONCOCT.out.versions.first() )
    }

    emit:
    metabat2depths   = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    maxbin_bins      = maxbin_output
    concoct_bins     = concoct_output
    metabat_bins     = metabat_output
    versions         = ch_versions
}