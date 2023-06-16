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
include { GZIP } from '../modules/utils'


process CHECKM_TABLE_FOR_DREP_GENOMES {
    tag "${name}"

    publishDir(
        path: "${params.outdir}/",
        mode: "copy"
    )

    input:
    tuple val(name), path(checkm), path(genomes_list)

    output:
    tuple val(name), path("checkm_results_MAGs.tab"), emit: checkm_results_MAGs

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
        cleaned_bins = CLEAN_AND_FILTER_BINS.out.filtered_bins.groupTuple()
        cleaned_bins.view()
        CHECKM2(cleaned_bins, ref_checkm.first())

        drep_input = cleaned_bins.combine(CHECKM2.out.checkm_table, by: 0)
        prok_drep_args = channel.value('-pa 0.9 -sa 0.95 -nc 0.6 -cm larger -comp 50 -con 5')
        DREP(drep_input, prok_drep_args, channel.value('prok'))

        // input: tuple(name, genomes/*, metabat_depth)
        metabat_depth = input_data.map(item -> tuple(item[0], item[2]))
        COVERAGE_RECYCLER(DREP.out.dereplicated_genomes.combine(metabat_depth, by: 0))

        CHANGE_UNDERSCORE_TO_DOT(DREP.out.dereplicated_genomes.transpose())

        DETECT_RRNA(DREP.out.dereplicated_genomes.transpose(), ref_rfam_rrna_models.first())

        GTDBTK(CHANGE_UNDERSCORE_TO_DOT.out.return_contigs.groupTuple(), ref_gtdbtk.first())

        CHECKM_TABLE_FOR_DREP_GENOMES(CHECKM2.out.checkm_table.combine(DREP.out.dereplicated_genomes_list, by:0))

        //GZIP(DREP.out.dereplicated_genomes.map(item -> item[1]))
}
