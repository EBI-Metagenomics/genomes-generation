/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Prokaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLEAN_AND_FILTER_BINS } from '../subworkflows/subwf_clean_and_filter_bins'
include { CHECKM2 } from '../modules/checkm2'
include { DREP } from '../modules/drep'
include { COVERAGE_RECYCLER } from '../modules/cov_recycler'
include { DETECT_RRNA } from '../modules/detect_rrna'
include { GTDBTK } from '../modules/gtdbtk'
include { CHANGE_UNDERSCORE_TO_DOT } from '../modules/utils'


process CHECKM_TABLE_FOR_DREP_GENOMES {

    publishDir(
        path: "${params.outdir}/",
        mode: "copy"
    )

    input:
    path(checkm)
    path(genomes_list)

    output:
    path("checkm_results_MAGs.tab"), emit: checkm_results_MAGs

    script:
    """
    grep -f ${genomes_list} ${checkm} > checkm_results_MAGs.tab
    """
}


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
        bins = CLEAN_AND_FILTER_BINS.out.filtered_bins.collect()

        CHECKM2(channel.value("aggregated"), bins, ref_checkm)

        prok_drep_args = channel.value('-pa 0.9 -sa 0.95 -nc 0.6 -cm larger -comp 50 -con 5')
        DREP(CHECKM2.out.checkm2_results, prok_drep_args, channel.value('prok'))

        metabat_depth = input_data.map(item -> item[2]).collectFile(name:"aggregated_metabat_depth.txt", skip:1, keepHeader:true)
        metabat_depth.view()
        COVERAGE_RECYCLER(DREP.out.dereplicated_genomes, metabat_depth)

        CHANGE_UNDERSCORE_TO_DOT(DREP.out.dereplicated_genomes.transpose())

        DETECT_RRNA(DREP.out.dereplicated_genomes.transpose(), ref_rfam_rrna_models.first())

        GTDBTK(CHANGE_UNDERSCORE_TO_DOT.out.return_files.groupTuple(), ref_gtdbtk.first())

        CHECKM_TABLE_FOR_DREP_GENOMES(CHECKM2.out.checkm2_results.map(item -> item[2]), DREP.out.dereplicated_genomes_list.map(item -> item[1]))

    emit:
        prok_mags = CHANGE_UNDERSCORE_TO_DOT.out.return_files.map(item -> item[1])
}
