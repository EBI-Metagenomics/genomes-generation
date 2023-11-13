/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Prokaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REFINEMENT as BIN_REFINEMENT       } from './mgbinrefinder/binrefinder'
include { CLEAN_AND_FILTER_BINS              } from './clean_and_filter_bins'

include { CHECKM2                            } from '../../modules/local/checkm2/main'
include { DREP                               } from '../../modules/local/drep/main'
include { COVERAGE_RECYCLER                  } from '../../modules/local/coverage_recycler/main'
include { DETECT_RRNA                        } from '../../modules/local/detect_rrna/main'
include { GTDBTK                             } from '../../modules/local/gtdbtk/main'
include { CHANGE_UNDERSCORE_TO_DOT           } from '../../modules/local/utils'


process CHECKM2_TABLE_FOR_DREP_GENOMES {

    input:
    path(checkm_filtered_genomes_dir)
    path(dereplicated_genomes_tsv)

    output:
    path("checkm_results_mags.tab"), emit: checkm_results_mags

    script:
    """
    grep -f ${dereplicated_genomes_tsv} ${checkm_filtered_genomes_dir} > checkm_results_mags.tab
    """
}


workflow PROK_MAGS_GENERATION {

    take:
    collected_binners_and_depth // tuple( meta, concoct, metabat, maxbin, depth_file)
    cat_db_folder
    cat_diamond_db
    cat_taxonomy_db
    gunc_db
    checkm2_db
    gtdbtk_db
    rfam_rrna_models

    main:

    ch_versions = Channel.empty()

    collected_binners = collected_binners_and_depth.map { meta, concot_bins, maxbin_bins, metabat_bins, _ -> 
        [ meta, concot_bins, maxbin_bins, metabat_bins ]
    }

    metabat_depth = collected_binners_and_depth.map { it -> it[4] }

    // -- bin refinement //
    BIN_REFINEMENT( collected_binners, checkm2_db )

    ch_versions.mix( BIN_REFINEMENT.out.versions.first() )

    // -- clean bins
    CLEAN_AND_FILTER_BINS( 
        BIN_REFINEMENT.out.bin_ref_bins,
        cat_db_folder,
        cat_diamond_db,
        cat_taxonomy_db,
        gunc_db
    )

    // -- aggregate bins by samples
    // -- checkm2 on ALL bins in all samples
    all_bins = CLEAN_AND_FILTER_BINS.out.bins.collect().map { it ->
        def meta = [:]
        meta.id = "aggregated"
        return tuple( meta, it )
    }

    CHECKM2( channel.value("aggregated"), all_bins, checkm2_db )

    // -- drep
    prok_drep_args = channel.value('-pa 0.9 -sa 0.95 -nc 0.6 -cm larger -comp 50 -con 5')

    DREP( CHECKM2.out.stats, prok_drep_args, channel.value('prokaryotes') )

    // -- coverage -- //
    COVERAGE_RECYCLER( DREP.out.dereplicated_genomes, metabat_depth.collect() )

    dereplicated_genomes = DREP.out.dereplicated_genomes.map { it -> it[1] }.flatten()

    CHANGE_UNDERSCORE_TO_DOT( dereplicated_genomes )

    // -- RNA -- //
    DETECT_RRNA( 
        dereplicated_genomes,
        rfam_rrna_models
    )

    // -- Taxonomy --//
    GTDBTK( CHANGE_UNDERSCORE_TO_DOT.out.return_files.collect(), gtdbtk_db )

    // -- checkm_results_mags.txt -- //
    // Both channels will have only one element
    CHECKM2_TABLE_FOR_DREP_GENOMES(
        CHECKM2.out.stats.map { map, bins, stats -> stats },
        DREP.out.dereplicated_genomes_list.map { meta, genomes_list_tsv -> genomes_list_tsv }
    )

    emit:
    prok_mags = CHANGE_UNDERSCORE_TO_DOT.out.return_files
    versions = ch_versions
}
