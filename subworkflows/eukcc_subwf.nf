include { EUKCC } from '../modules/eukcc'
include { LINKTABLE } from '../modules/eukcc'
/*
    ~~~~~~~~~~~~~~~~~~~~~~
     Run EukCC subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~
*/
workflow EUKCC_SUBWF {
    take:
        name
        bins
        bam
        bam_index
        eukcc_db
    main:
        LINKTABLE(name, bins, bam, bam_index)
        EUKCC(name, LINKTABLE.out.links_table, eukcc_db, bins)
    emit:
        eukcc_result = EUKCC.out.eukcc_results
}
