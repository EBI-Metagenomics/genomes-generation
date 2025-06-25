include { GUNZIP as GUNZIP_ASSEMBLY   } from '../../modules/local/utils'

include { CONCOCT_CONCOCT              } from '../../modules/nf-core/concoct/concoct'
include { CONCOCT_MERGECUTUPCLUSTERING } from '../../modules/nf-core/concoct/mergecutupclustering'
include { CONCOCT_EXTRACTFASTABINS     } from '../../modules/nf-core/concoct/extractfastabins'


workflow CONCOCT_SUBWF {

    take:
    ch_fasta_tsv // channel (mandatory): [ val(meta), tsv, concoct_fasta, assembly_fasta ]

    main:
    ch_versions = Channel.empty()

    /*
    * --- uncompress assembly fasta ---
    */
    GUNZIP_ASSEMBLY(
        ch_fasta_tsv.map { meta, _1, _2, assembly -> [meta, assembly] }
    )

    CONCOCT_CONCOCT( 
        ch_fasta_tsv.map{ meta, tsv, concoct_fasta, assembly_fasta -> [meta, tsv, concoct_fasta] } 
    )

    ch_versions = ch_versions.mix(CONCOCT_CONCOCT.out.versions.first())

    CONCOCT_MERGECUTUPCLUSTERING ( 
        CONCOCT_CONCOCT.out.clustering_csv 
    )

    ch_versions = ch_versions.mix( CONCOCT_MERGECUTUPCLUSTERING.out.versions.first())

    ch_mergecutupclustering_for_extractfastabins = 

    CONCOCT_EXTRACTFASTABINS ( 
        GUNZIP_ASSEMBLY.out.uncompressed
        .join(CONCOCT_MERGECUTUPCLUSTERING.out.csv, failOnMismatch: true)
    )

    ch_versions = ch_versions.mix(CONCOCT_EXTRACTFASTABINS.out.versions.first())

    emit:

    bins     = CONCOCT_EXTRACTFASTABINS.out.fasta       // channel: [ val(meta), [ fasta ] ]
    versions = ch_versions                                         // channel: [ versions.yml ]
}