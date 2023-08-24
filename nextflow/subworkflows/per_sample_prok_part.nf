/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Prokaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REFINEMENT as BIN_REFINEMENT       } from './mgbinrefinder/binrefinder'
include { CLEAN_AND_FILTER_BINS              } from '../subworkflows/local/subwf_clean_and_filter_bins'
include { CHECKM2                            } from '../modules/local/checkm2/main'
include { DREP                               } from '../modules/local/drep/main'
include { COVERAGE_RECYCLER                  } from '../modules/local/coverage_recycler/main'
include { DETECT_RRNA                        } from '../modules/local/detect_rrna/main'
include { GTDBTK                             } from '../modules/local/gtdbtk/main'
include { CHANGE_UNDERSCORE_TO_DOT           } from '../modules/local/utils'


process CHECKM_TABLE_FOR_DREP_GENOMES {

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
        collected_binners_and_depth
        ref_catdb
        ref_cat_diamond
        ref_cat_taxonomy
        ref_gunc
        ref_checkm
        ref_gtdbtk
        ref_rfam_rrna_models
    main:
        collected_binners = collected_binners_and_depth.map{it -> [it[0], it[1], it[2], it[3]]}
        metabat_depth = collected_binners_and_depth.map{it -> it[4]}

        // -- bin refinement
        BIN_REFINEMENT(collected_binners, ref_checkm)
        BIN_REFINEMENT.out.bin_ref_bins.view()
        // -- clean bins
        CLEAN_AND_FILTER_BINS(BIN_REFINEMENT.out.bin_ref_bins, ref_catdb.first(), ref_cat_diamond.first(), ref_cat_taxonomy.first(), ref_gunc.first())

        // -- aggregate bins by samples
        // -- checkm2 on ALL bins in all samples
        all_bins = CLEAN_AND_FILTER_BINS.out.bins.collect().map{ it ->
                                                                    def meta = [:]
                                                                    meta.id = "aggregated"
                                                                    return tuple(meta, it)}
        CHECKM2(channel.value("aggregated"), all_bins, ref_checkm)

        // -- drep
        prok_drep_args = channel.value('-pa 0.9 -sa 0.95 -nc 0.6 -cm larger -comp 50 -con 5')
        DREP(CHECKM2.out.checkm2_results, prok_drep_args, channel.value('prok'))

        // -- coverage
        COVERAGE_RECYCLER(DREP.out.dereplicated_genomes, metabat_depth.collect())

        CHANGE_UNDERSCORE_TO_DOT(DREP.out.dereplicated_genomes.map{it->it[1]}.flatten())

        // -- RNA
        DETECT_RRNA(DREP.out.dereplicated_genomes.map{it->it[1]}.flatten(), ref_rfam_rrna_models.first())

        // -- Taxonomy
        GTDBTK(CHANGE_UNDERSCORE_TO_DOT.out.return_files.collect(), ref_gtdbtk)

        // -- checkm_results_MAGs.txt
        CHECKM_TABLE_FOR_DREP_GENOMES(CHECKM2.out.checkm2_results.map(item -> item[2]), DREP.out.dereplicated_genomes_list.map(item -> item[1]))

    emit:
        prok_mags = CHANGE_UNDERSCORE_TO_DOT.out.return_files
}
