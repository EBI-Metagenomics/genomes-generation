/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Prokaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REFINEMENT as BIN_REFINEMENT       } from './mgbinrefinder/binrefinder'
include { CLEAN_AND_FILTER_BINS              } from '../subworkflows/local/clean_and_filter_bins'

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


workflow PROK_MAGS_GENERATION {
    take:
    collected_binners_and_depth // //
    cat_db
    cat_diamond_db
    cat_taxonomy_db
    gunc_db
    checkm2_db
    gtdbtk_db
    rfam_rrna_models

    main:
    collected_binners = collected_binners_and_depth.map{ it -> [it[0], it[1], it[2], it[3]] }
    metabat_depth = collected_binners_and_depth.map{ it -> it[4] }

    // -- bin refinement //
    BIN_REFINEMENT( collected_binners, checkm2_db )

    // -- clean bins
    CLEAN_AND_FILTER_BINS( BIN_REFINEMENT.out.bin_ref_bins, cat_db, cat_diamond_db, cat_taxonomy_db, gunc_db )

    // -- aggregate bins by samples
    // -- checkm2 on ALL bins in all samples
    all_bins = CLEAN_AND_FILTER_BINS.out.bins.collect().map{ it ->
                                                                def meta = [:]
                                                                meta.id = "aggregated"
                                                                return tuple(meta, it)
                                                            }
    CHECKM2( channel.value("aggregated"), all_bins, checkm2_db )

    // -- drep
    prok_drep_args = channel.value('-pa 0.9 -sa 0.95 -nc 0.6 -cm larger -comp 50 -con 5')

    DREP( CHECKM2.out.checkm2_results, prok_drep_args, channel.value('prok') )

    // -- coverage -- //
    COVERAGE_RECYCLER( DREP.out.dereplicated_genomes, metabat_depth.collect() )

    CHANGE_UNDERSCORE_TO_DOT(DREP.out.dereplicated_genomes.map{ it -> it[1] }.flatten())

    // -- RNA -- //
    DETECT_RRNA( 
        DREP.out.dereplicated_genomes.map{ it -> it[1] }.flatten(),
        rfam_rrna_models
    )

    // -- Taxonomy --//
    GTDBTK( CHANGE_UNDERSCORE_TO_DOT.out.return_files.collect(), gtdbtk_db )

    // -- checkm_results_MAGs.txt -- //
    CHECKM_TABLE_FOR_DREP_GENOMES(
        CHECKM2.out.checkm2_results.map { it -> it[2] },
        DREP.out.dereplicated_genomes_list.map { it -> it[1] }
    )

    emit:
    prok_mags = CHANGE_UNDERSCORE_TO_DOT.out.return_files
}
