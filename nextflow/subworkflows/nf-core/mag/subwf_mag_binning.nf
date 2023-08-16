/*
 * Binning with MetaBAT2, MaxBin2 and Concoct
 */
include { CREATE_FOLDER as FOLDER_METABAT       } from '../../../modules/local/utils'
include { CREATE_FOLDER as FOLDER_CONCOCT       } from '../../../modules/local/utils'
include { CREATE_FOLDER as FOLDER_MAXBIN        } from '../../../modules/local/utils'

include { METABAT2_METABAT2                     } from '../../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS  } from '../../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { MAXBIN2                               } from '../../../modules/nf-core/maxbin2/main'

include { CONVERT_DEPTHS                        } from '../../../modules/nf-core/mag/convert_depths'
include { ADJUST_MAXBIN2_EXT                    } from '../../../modules/nf-core/mag/adjust_maxbin2_ext'
include { FASTA_BINNING_CONCOCT                 } from '../fasta_binning_concoct/main'

/*
 * Get number of columns in file (first line)
 */
def getColNo(filename) {
    lines  = file(filename).readLines()
    return lines[0].split('\t').size()
}

workflow BINNING {
    take:
    assemblies           // channel: [ val(meta), path(assembly), path(bams), path(bais) ]
    reads                // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()

    // generate coverage depths for each contig
    ch_summarizedepth_input = assemblies
                                .map { meta, assembly, bams, bais ->
                                        def meta_new = meta.clone()
                                    [ meta_new, bams, bais ]
                                }

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )

    ch_metabat_depths = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, depths ]
        }

    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first())

    // combine depths back with assemblies
    ch_metabat2_input = assemblies
        .map { meta, assembly, bams, bais ->
            def meta_new = meta.clone()
            meta_new['binner'] = 'MetaBAT2'

            [ meta_new, assembly, bams, bais ]
        }
        .join( ch_metabat_depths, by: 0 )
        .map { meta, assembly, bams, bais, depths ->
            [ meta, assembly, depths ]
        }

    // convert metabat2 depth files to maxbin2
    if ( !params.skip_maxbin2 ) {
        CONVERT_DEPTHS ( ch_metabat2_input )
        ch_maxbin2_input = CONVERT_DEPTHS.out.output
            .map { meta, assembly, reads, depth ->
                    def meta_new = meta.clone()
                    meta_new['binner'] = 'MaxBin2'

                [ meta_new, assembly, reads, depth ]
            }
        ch_versions = ch_versions.mix(CONVERT_DEPTHS.out.versions.first())
    }

    // final gzipped bins
    // run binning
    if ( !params.skip_metabat2 ) {
        METABAT2_METABAT2 ( ch_metabat2_input )
        // before decompressing first have to separate and re-group due to limitation of GUNZIP module
        size_binner = METABAT2_METABAT2.out.fasta.map{it -> [it[0], it[1].collect().size()]}
        size_binner.view()
        ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions.first())
    }
    if ( !params.skip_maxbin2 ) {
        MAXBIN2 ( ch_maxbin2_input )
        ADJUST_MAXBIN2_EXT ( MAXBIN2.out.binned_fastas )
        size_binner = ADJUST_MAXBIN2_EXT.out.renamed_bins.map{it -> [it[0], it[1].collect().size()]}
        size_binner.view()
        ch_versions = ch_versions.mix(MAXBIN2.out.versions)
    }
    if ( !params.skip_concoct ){

        ch_concoct_input = assemblies
                            .map { meta, bins, bams, bais ->
                                def meta_new = meta.clone()
                                meta_new['binner'] = 'CONCOCT'
                                [ meta_new, bins, bams, bais ]
                            }
                            .multiMap {
                                meta, bins, bams, bais ->
                                    bins: [ meta, bins ]
                                    bams: [ meta, bams, bais ]
                            }

        FASTA_BINNING_CONCOCT ( ch_concoct_input )
        size_binner = FASTA_BINNING_CONCOCT.out.bins.map{it -> [it[0], it[1].collect().size()]}
        size_binner.view()
        ch_versions = ch_versions.mix(FASTA_BINNING_CONCOCT.out.versions)
    }

    assembly = assemblies.map{it -> [it[0],it[1]]}
    meta = assemblies.map{it -> it[0]}
    FOLDER_MAXBIN(ADJUST_MAXBIN2_EXT.out.renamed_bins.map{meta, bins ->
                                                            meta.remove('binner')
                                                            [meta, bins]
                                                         }, channel.value('maxbin'))
    FOLDER_CONCOCT(FASTA_BINNING_CONCOCT.out.bins.map{meta, bins ->
                                                            meta.remove('binner')
                                                            [meta, bins]
                                                         }, channel.value('concoct'))
    FOLDER_METABAT(METABAT2_METABAT2.out.fasta.map{meta, bins ->
                                                            meta.remove('binner')
                                                            [meta, bins]
                                                         }, channel.value('metabat'))
    emit:
        concoct_bins                                 = FOLDER_CONCOCT.out.output
        metabat_bins                                 = FOLDER_METABAT.out.output
        maxbin_bins                                  = FOLDER_MAXBIN.out.output
        metabat2depths                               = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        versions                                     = ch_versions
}