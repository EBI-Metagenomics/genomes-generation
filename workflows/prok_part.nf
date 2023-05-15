/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Prokaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLEAN_AND_FILTER_BINS } from '../subworkflows/clean_and_filter_bins'
include { CHECKM2 } from '../modules/checkm2'
include { DREP } from '../modules/drep'

workflow PROK_SUBWF {
    take:
        sample_name
        bins
        metabat_depth
        ref_catdb
        ref_cat_diamond
        ref_cat_taxonomy
        ref_gunc
        ref_checkm
        ref_gtdbtk
        ref_rfam_rrna_models
    main:
        CLEAN_AND_FILTER_BINS(sample_name, bins, ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc)

        CHECKM2(CLEAN_AND_FILTER_BINS.out.filtered_bins, ref_checkm)

        nc_drep = channel.value(0.60)
        DREP(CLEAN_AND_FILTER_BINS.out.filtered_bins, CHECKM2.out.checkm_table, nc_drep)

        DETECT_RRNA(DREP.out.dereplicated_genomes, ref_rfam_rrna_models)

        COVERAGE_RECYCLER(DREP.out.dereplicated_genomes, metabat_depth)

        CHANGE_UNDERSCORE_TO_DOT(DREP.out.dereplicated_genomes)
        GTDBTK(CHANGE_UNDERSCORE_TO_DOT.out.return_contigs, ref_gtdbtk)
    emit:
        genomes = CHANGE_UNDERSCORE_TO_DOT.out.return_contigs
}
