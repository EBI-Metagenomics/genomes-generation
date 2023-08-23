include { CHECKM2 as CHECKM_1                                               } from '../../modules/local/checkm2/main'
include { CHECKM2 as CHECKM_2                                               } from '../../modules/local/checkm2/main'
include { CHECKM2 as CHECKM_3                                               } from '../../modules/local/checkm2/main'
include { CHECKM2 as CHECKM_FINAL                                           } from '../../modules/local/checkm2/main'
include { RENAME_AND_CHECK_SIZE_BINS as RENAME_AND_CHECK_SIZE_BINS_BINNER1  } from '../../modules/mgbinrefinder/utils'
include { RENAME_AND_CHECK_SIZE_BINS as RENAME_AND_CHECK_SIZE_BINS_BINNER2  } from '../../modules/mgbinrefinder/utils'
include { RENAME_AND_CHECK_SIZE_BINS as RENAME_AND_CHECK_SIZE_BINS_BINNER3  } from '../../modules/mgbinrefinder/utils'
include { REFINE as REFINE12                                                } from './refine'
include { REFINE as REFINE13                                                } from './refine'
include { REFINE as REFINE23                                                } from './refine'
include { REFINE as REFINE123                                               } from './refine'
include { CONSOLIDATE_BINS                                                  } from '../../modules/mgbinrefinder/consolidate_bins'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/


workflow REFINEMENT {
    take:
        collected_binners
        ref_checkm
    main:
        meta = collected_binners.map{it -> it[0]}
        binner1 = collected_binners.map{it -> [it[0], it[1]]}.transpose()
        binner2 = collected_binners.map{it -> [it[0], it[2]]}.transpose()
        binner3 = collected_binners.map{it -> [it[0], it[3]]}.transpose()

        RENAME_AND_CHECK_SIZE_BINS_BINNER1(channel.value("binner1"), binner1)
        RENAME_AND_CHECK_SIZE_BINS_BINNER2(channel.value("binner2"), binner2)
        RENAME_AND_CHECK_SIZE_BINS_BINNER3(channel.value("binner3"), binner3)

        // collect by meta
        renamed_binner1 = RENAME_AND_CHECK_SIZE_BINS_BINNER1.out.renamed.groupTuple(by:0)
        renamed_binner2 = RENAME_AND_CHECK_SIZE_BINS_BINNER2.out.renamed.groupTuple(by:0)
        renamed_binner3 = RENAME_AND_CHECK_SIZE_BINS_BINNER3.out.renamed.groupTuple(by:0)

        size_1 = renamed_binner1.map{it -> it[1]}.collect().size()
        //size_1 = renamed_binner1.size()
        //size_1.subscribe { println "Filter results: Binner1: $it.value" }
        //renamed_binner2 = collected_binner2.renamed.collect()
        //size_2 = renamed_binner2.size()
        //size_2.subscribe { println "Filter results: Binner2: $it.value" }
        //renamed_binner3 = collected_binner3.renamed.collect()
        //size_3 = renamed_binner3.size()
        //size_3.subscribe { println "Filter results: Binner3: $it.value" }

        REFINE12(channel.value("binner12"), renamed_binner1, renamed_binner2, false, ref_checkm)
        REFINE13(channel.value("binner13"), renamed_binner1, renamed_binner3, false, ref_checkm)
        REFINE23(channel.value("binner23"), renamed_binner2, renamed_binner3, false, ref_checkm)
        REFINE123(channel.value("binner123"), renamed_binner1, renamed_binner2, renamed_binner3, ref_checkm)

        CHECKM_1(channel.value("binner1"), renamed_binner1, ref_checkm)
        CHECKM_2(channel.value("binner2"), renamed_binner2, ref_checkm)
        CHECKM_3(channel.value("binner3"), renamed_binner3, ref_checkm)

        binners = CHECKM_1.out.checkm2_results_filtered.combine(
                    CHECKM_2.out.checkm2_results_filtered, by:0).combine(
                    CHECKM_3.out.checkm2_results_filtered,by:0).combine(
                    REFINE12.out.filtered_bins,by:0).combine(
                    REFINE13.out.filtered_bins,by:0).combine(
                    REFINE23.out.filtered_bins,by:0).combine(
                    REFINE123.out.filtered_bins,by:0)

        stats = CHECKM_1.out.checkm2_results_filtered_stats.combine(
                    CHECKM_2.out.checkm2_results_filtered_stats, by:0).combine(
                    CHECKM_3.out.checkm2_results_filtered_stats, by:0).combine(
                    REFINE12.out.filtered_bins_stats, by:0).combine(
                    REFINE13.out.filtered_bins_stats, by:0).combine(
                    REFINE23.out.filtered_bins_stats, by:0).combine(
                    REFINE123.out.filtered_bins_stats, by:0)

        CONSOLIDATE_BINS(binners, stats)

        //CHECKM_FINAL(channel.value("final"), CONSOLIDATE_BINS.out.dereplicated_bins, ref_checkm)

    emit:
        //checkm_stats = CHECKM_FINAL.out.checkm_results
        bin_ref_bins = CONSOLIDATE_BINS.out.dereplicated_bins
}
