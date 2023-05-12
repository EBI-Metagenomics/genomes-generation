/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Prokaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLEAN_AND_FILTER_BINS } from '../subworkflows/clean_and_filter_bins'
include { CHECKM_SUBWF } from '../subworkflows/checkm_subwf'
include { DREP } from '../modules/drep'

workflow PROK_SUBWF {
    take:
        sample_name
        bins
        ref_catdb
        ref_cat_diamond
        ref_cat_taxonomy
        ref_gunc
        ref_checkm
    main:
        CLEAN_AND_FILTER_BINS(sample_name, bins, ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc)

        CHECKM_SUBWF(sample_name, CLEAN_AND_FILTER_BINS.out.filtered_bins, ref_checkm)

        nc_drep = channel.value(0.60)
        DREP(CLEAN_AND_FILTER_BINS.out.filtered_bins, CHECKM_SUBWF.out.checkm_table, nc_drep)
    emit:

}
