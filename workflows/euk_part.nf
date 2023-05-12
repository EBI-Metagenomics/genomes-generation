/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EUKCC_SUBWF } from '../subworkflows/eukcc_subwf'

workflow EUK_SUBWF {
    take:
        mode
        sample_name
        bins

    main:
        EUKCC_SUBWF()

    emit:

}