/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Prokaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLEAN_AND_FILTER_BINS } from '../subworkflows/subwf_clean_and_filter_bins'
include { CHECKM2 } from '../modules/checkm2'
include { DREP } from '../modules/drep'
include { COVERAGE_RECYCLER } from '../modules/cov_recycler'
include { CHANGE_UNDERSCORE_TO_DOT } from '../modules/utils'

workflow PROK_SUBWF {
    take:
        input_data                // tuple(accession, [bins .fa, ...], file_depth_metabat2)
        ref_catdb
        ref_cat_diamond
        ref_cat_taxonomy
        ref_gunc
        ref_checkm
        ref_gtdbtk
        ref_rfam_rrna_models
    main:
        bins = input_data.map(item -> tuple(item[0], item[1]))
        CLEAN_AND_FILTER_BINS(bins, ref_catdb.first(), ref_cat_diamond.first(), ref_cat_taxonomy.first(), ref_gunc.first())

        CHECKM2(CLEAN_AND_FILTER_BINS.out.filtered_bins, ref_checkm.first())

        drep_input = CLEAN_AND_FILTER_BINS.out.filtered_bins.combine(CHECKM2.out.checkm_table, by: 0)
        prok_drep_args = channel.value('-pa 0.9 -sa 0.95 -nc 0.6 -cm larger -comp 50 -con 5')
        DREP(drep_input, prok_drep_args)

        // input: tuple(name, genomes/*, metabat_depth)
        metabat_depth = input_data.map(item -> tuple(item[0], item[2]))
        COVERAGE_RECYCLER(DREP.out.dereplicated_genomes.combine(metabat_depth, by: 0))

        CHANGE_UNDERSCORE_TO_DOT(DREP.out.dereplicated_genomes.transpose())

        //DETECT_RRNA(DREP.out.dereplicated_genomes, ref_rfam_rrna_models)

        //GTDBTK(CHANGE_UNDERSCORE_TO_DOT.out.return_contigs, ref_gtdbtk)
    emit:
        out = DREP.out.dereplicated_genomes
        //genomes = CHANGE_UNDERSCORE_TO_DOT.out.return_contigs
}
