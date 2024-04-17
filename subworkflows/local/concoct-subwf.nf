include { CONCOCT_CONCOCT              } from '../../modules/nf-core/concoct/concoct/main.nf'
include { CONCOCT_MERGECUTUPCLUSTERING } from '../../modules/nf-core/concoct/mergecutupclustering/main.nf'
include { CONCOCT_EXTRACTFASTABINS     } from '../../modules/nf-core/concoct/extractfastabins/main.nf'

workflow CONCOCT_SUBWF {

    take:
    ch_fasta_tsv // channel (mandatory): [ val(meta), tsv, fasta ]

    main:
    ch_versions = Channel.empty()

    ch_fasta = ch_fasta_tsv.map{ meta, fasta, tsv -> [meta, fasta] }

    CONCOCT_CONCOCT( ch_fasta_tsv )

    ch_versions = ch_versions.mix(CONCOCT_CONCOCT.out.versions.first())

    CONCOCT_MERGECUTUPCLUSTERING ( CONCOCT_CONCOCT.out.clustering_csv )

    ch_versions = ch_versions.mix( CONCOCT_MERGECUTUPCLUSTERING.out.versions.first())

    ch_mergecutupclustering_for_extractfastabins = ch_fasta
                                                    .join(CONCOCT_MERGECUTUPCLUSTERING.out.csv, failOnMismatch: true)

    CONCOCT_EXTRACTFASTABINS ( ch_mergecutupclustering_for_extractfastabins )

    ch_versions = ch_versions.mix(CONCOCT_EXTRACTFASTABINS.out.versions.first())

    emit:

    bins                = CONCOCT_EXTRACTFASTABINS.out.fasta       // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions                                         // channel: [ versions.yml ]
}