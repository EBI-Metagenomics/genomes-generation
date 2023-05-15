/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EUKCC } from '../modules/eukcc'
include { LINKTABLE } from '../modules/eukcc'

workflow EUK_SUBWF {
    take:
        mode
        name
        bins
        bam
        bam_index

    main:
        LINKTABLE(name, bins, bam, bam_index)
        EUKCC(name, LINKTABLE.out.links_table, eukcc_db, bins)

    emit:

}