include { CHECKM2 as CHECKM_1                           } from '../../modules/checkm2'
include { CHECKM2 as CHECKM_2                           } from '../../modules/checkm2'
include { CHECKM2 as CHECKM_3                           } from '../../modules/checkm2'
include { CHECKM2 as CHECKM_FINAL                       } from '../../modules/checkm2'
include { RENAME_AND_CHECK_SIZE_BINS as RENAME_BINNER1  } from '../..modules/mgbinrefinder/utils'
include { RENAME_AND_CHECK_SIZE_BINS as RENAME_BINNER2  } from '../../modules/mgbinrefinder/utils'
include { RENAME_AND_CHECK_SIZE_BINS as RENAME_BINNER3  } from '../../modules/mgbinrefinder/utils'
include { REFINE as REFINE12                            } from './refine'
include { REFINE as REFINE13                            } from './refine'
include { REFINE as REFINE23                            } from './refine'
include { REFINE as REFINE123                           } from './refine'
include { CONSOLIDATE_BINS                              } from '../../modules/mgbinrefinder/consolidate_bins'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow REFINEMENT {
    RENAME_BINNER1(channel.value("binner1"), binner1)
    RENAME_BINNER2(channel.value("binner2"), binner2)
    RENAME_BINNER3(channel.value("binner3"), binner3)

    renamed_binner1 = RENAME_BINNER1.out.renamed.collect()
    size_1 = renamed_binner1.size()
    size_1.subscribe { println "Filter results: Binner1: $it.value" }
    renamed_binner2 = RENAME_BINNER2.out.renamed.collect()
    size_2 = renamed_binner2.size()
    size_2.subscribe { println "Filter results: Binner2: $it.value" }
    renamed_binner3 = RENAME_BINNER3.out.renamed.collect()
    size_3 = renamed_binner3.size()
    size_3.subscribe { println "Filter results: Binner3: $it.value" }

    REFINE12(channel.value("binner12"), renamed_binner1, renamed_binner2, channel.fromPath('NO_FILE'), ref_checkm)
    REFINE13(channel.value("binner13"), renamed_binner1, renamed_binner3, channel.fromPath('NO_FILE'), ref_checkm)
    REFINE23(channel.value("binner23"), renamed_binner2, renamed_binner3, channel.fromPath('NO_FILE'), ref_checkm)
    REFINE123(channel.value("binner123"), renamed_binner1, renamed_binner2, renamed_binner3, ref_checkm)

    CHECKM_1(channel.value("binner1"), renamed_binner1, ref_checkm)
    CHECKM_2(channel.value("binner2"), renamed_binner2, ref_checkm)
    CHECKM_3(channel.value("binner3"), renamed_binner3, ref_checkm)

    binners = CHECKM_1.out.checkm2_results_filtered.concat(
                CHECKM_2.out.checkm2_results_filtered).concat(
                CHECKM_3.out.checkm2_results_filtered).concat(
                REFINE12.out.filtered_bins).concat(
                REFINE13.out.filtered_bins).concat(
                REFINE23.out.filtered_bins).concat(
                REFINE123.out.filtered_bins)
    stats = CHECKM_1.out.checkm2_results_filtered_stats.concat(
                CHECKM_2.out.checkm2_results_filtered_stats).concat(
                CHECKM_3.out.checkm2_results_filtered_stats).concat(
                REFINE12.out.filtered_bins_stats).concat(
                REFINE13.out.filtered_bins_stats).concat(
                REFINE23.out.filtered_bins_stats).concat(
                REFINE123.out.filtered_bins_stats)
    CONSOLIDATE_BINS(binners.collect(), stats.collect())

    CHECKM_FINAL(channel.value("final"), CONSOLIDATE_BINS.out.dereplicated_bins, ref_checkm)
}
