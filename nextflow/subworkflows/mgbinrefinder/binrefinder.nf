include { CHECKM2 as CHECKM2_BINNER1                                        } from '../../modules/local/checkm2/main'
include { CHECKM2 as CHECKM2_BINNER2                                        } from '../../modules/local/checkm2/main'
include { CHECKM2 as CHECKM2_BINNER3                                        } from '../../modules/local/checkm2/main'
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
        collected_binners  // [meta, concoct, maxbin, metabat]
        ref_checkm
    main:
        meta = collected_binners.map{it -> it[0]}

        // if it needs to be a list of (meta, bin1), (meta, bin2) - use .transpose()
        binner1 = collected_binners.map{it -> [it[0], it[1]]}
        binner2 = collected_binners.map{it -> [it[0], it[2]]}
        binner3 = collected_binners.map{it -> [it[0], it[3]]}

        // TODO: it should be possible to filter by size using findAll() ot filter()
        // but I don't know how to do it properly with list of paths

        RENAME_AND_CHECK_SIZE_BINS_BINNER1(channel.value("binner1"), binner1)
        RENAME_AND_CHECK_SIZE_BINS_BINNER2(channel.value("binner2"), binner2)
        RENAME_AND_CHECK_SIZE_BINS_BINNER3(channel.value("binner3"), binner3)

        // collect by meta
        renamed_binner1 = RENAME_AND_CHECK_SIZE_BINS_BINNER1.out.renamed
        renamed_binner2 = RENAME_AND_CHECK_SIZE_BINS_BINNER2.out.renamed
        renamed_binner3 = RENAME_AND_CHECK_SIZE_BINS_BINNER3.out.renamed

        REFINE12(channel.value("binner12"), renamed_binner1, renamed_binner2, false, ref_checkm)
        REFINE13(channel.value("binner13"), renamed_binner1, renamed_binner3, false, ref_checkm)
        REFINE23(channel.value("binner23"), renamed_binner2, renamed_binner3, false, ref_checkm)
        REFINE123(channel.value("binner123"), renamed_binner1, renamed_binner2, renamed_binner3, ref_checkm)

        CHECKM2_BINNER1(channel.value("binner1"), renamed_binner1.map{it -> [it[0], it[1]]}, ref_checkm.first())
        CHECKM2_BINNER2(channel.value("binner2"), renamed_binner2.map{it -> [it[0], it[1]]}, ref_checkm.first())
        CHECKM2_BINNER3(channel.value("binner3"), renamed_binner3.map{it -> [it[0], it[1]]}, ref_checkm.first())

        binners = CHECKM2_BINNER1.out.checkm2_results_filtered.combine(
                    CHECKM2_BINNER2.out.checkm2_results_filtered, by:0).combine(
                    CHECKM2_BINNER3.out.checkm2_results_filtered,by:0).combine(
                    REFINE12.out.filtered_bins,by:0).combine(
                    REFINE13.out.filtered_bins,by:0).combine(
                    REFINE23.out.filtered_bins,by:0).combine(
                    REFINE123.out.filtered_bins,by:0)

        stats = CHECKM2_BINNER1.out.checkm2_results_filtered_stats.combine(
                    CHECKM2_BINNER2.out.checkm2_results_filtered_stats, by:0).combine(
                    CHECKM2_BINNER3.out.checkm2_results_filtered_stats, by:0).combine(
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
