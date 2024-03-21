include { CHECKM2 as CHECKM2_BINNER1                                        } from '../../../modules/local/checkm2/main'
include { CHECKM2 as CHECKM2_BINNER2                                        } from '../../../modules/local/checkm2/main'
include { CHECKM2 as CHECKM2_BINNER3                                        } from '../../../modules/local/checkm2/main'
include { RENAME_AND_CHECK_SIZE_BINS as RENAME_AND_CHECK_SIZE_BINS_BINNER1  } from '../../../modules/local/mgbinrefinder/utils'
include { RENAME_AND_CHECK_SIZE_BINS as RENAME_AND_CHECK_SIZE_BINS_BINNER2  } from '../../../modules/local/mgbinrefinder/utils'
include { RENAME_AND_CHECK_SIZE_BINS as RENAME_AND_CHECK_SIZE_BINS_BINNER3  } from '../../../modules/local/mgbinrefinder/utils'
include { CONSOLIDATE_BINS                                                  } from '../../../modules/local/mgbinrefinder/consolidate_bins'
include { REFINE as REFINE12                                                } from './refine'
include { REFINE as REFINE13                                                } from './refine'
include { REFINE as REFINE23                                                } from './refine'
include { REFINE as REFINE123                                               } from './refine'

workflow REFINEMENT {
    take:
    collected_binners
    ref_checkm

    main:

    ch_versions = Channel.empty()
    ch_log = Channel.empty()

    binner1 = collected_binners.map { meta, concot_bins, maxbin_bins, metabat_bins ->
        [ meta, concot_bins ]
    }

    binner2 = collected_binners.map { meta, concot_bins, maxbin_bins, metabat_bins -> 
        [ meta, maxbin_bins ]
    }

    binner3 = collected_binners.map { meta, concot_bins, maxbin_bins, metabat_bins ->
        [ meta, metabat_bins ]
    }

    RENAME_AND_CHECK_SIZE_BINS_BINNER1( "binner1", binner1 )
    RENAME_AND_CHECK_SIZE_BINS_BINNER2( "binner2", binner2 )
    RENAME_AND_CHECK_SIZE_BINS_BINNER3( "binner3", binner3 )

    // collect by meta
    renamed_binner1 = RENAME_AND_CHECK_SIZE_BINS_BINNER1.out.renamed
    renamed_binner2 = RENAME_AND_CHECK_SIZE_BINS_BINNER2.out.renamed
    renamed_binner3 = RENAME_AND_CHECK_SIZE_BINS_BINNER3.out.renamed

    binners_all = renamed_binner1.join(renamed_binner2, remainder: true).join(renamed_binner3, remainder: true)
    // replace null (generated with remainder to [])
    binners = binners_all.map{ meta, b1, b2, b3 ->
                        result = [meta]
                        if (b1) { result.add(b1) } else { result.add([]) }
                        if (b2) { result.add(b2) } else { result.add([]) }
                        if (b3) { result.add(b3) } else { result.add([]) }
                        return result
    }
    refine12_input = binners.map{meta, b1, b2, b3 -> [meta, b1, b2, []]}
    refine13_input = binners.map{meta, b1, b2, b3 -> [meta, b1, b3, []]}
    refine23_input = binners.map{meta, b1, b2, b3 -> [meta, b2, b3, []]}
    REFINE12( "binner12", refine12_input, ref_checkm )
    REFINE13( "binner13", refine13_input, ref_checkm )
    REFINE23( "binner23", refine23_input, ref_checkm )
    REFINE123( "binner123", binners, ref_checkm )

    CHECKM2_BINNER1( "binner1", renamed_binner1, ref_checkm )
    CHECKM2_BINNER2( "binner2", renamed_binner2, ref_checkm )
    CHECKM2_BINNER3( "binner3", renamed_binner3, ref_checkm )

    ch_versions = ch_versions.mix( CHECKM2_BINNER1.out.versions.first() )
    ch_versions = ch_versions.mix( CHECKM2_BINNER2.out.versions.first() )
    ch_versions = ch_versions.mix( CHECKM2_BINNER3.out.versions.first() )

    binners = CHECKM2_BINNER1.out.filtered_genomes
        .join(CHECKM2_BINNER2.out.filtered_genomes)
        .join(CHECKM2_BINNER3.out.filtered_genomes)
        .join(REFINE12.out.filtered_bins)
        .join(REFINE13.out.filtered_bins)
        .join(REFINE23.out.filtered_bins)
        .join(REFINE123.out.filtered_bins)

    stats = CHECKM2_BINNER1.out.filtered_stats
        .join(CHECKM2_BINNER2.out.filtered_stats)
        .join(CHECKM2_BINNER3.out.filtered_stats)
        .join(REFINE12.out.filtered_bins_stats)
        .join(REFINE13.out.filtered_bins_stats)
        .join(REFINE23.out.filtered_bins_stats)
        .join(REFINE123.out.filtered_bins_stats)

    consolidate_input = binners.join(stats)

    CONSOLIDATE_BINS( consolidate_input )
    ch_versions = ch_versions.mix( CONSOLIDATE_BINS.out.versions.first() )

    // collect logging
    ch_log = ch_log.mix( RENAME_AND_CHECK_SIZE_BINS_BINNER1.out.progress_log )
    ch_log = ch_log.mix( RENAME_AND_CHECK_SIZE_BINS_BINNER2.out.progress_log )
    ch_log = ch_log.mix( RENAME_AND_CHECK_SIZE_BINS_BINNER3.out.progress_log )
    ch_log = ch_log.mix( REFINE12.out.progress_log )
    ch_log = ch_log.mix( REFINE13.out.progress_log )
    ch_log = ch_log.mix( REFINE23.out.progress_log )
    ch_log = ch_log.mix( REFINE123.out.progress_log )
    ch_log = ch_log.mix( CHECKM2_BINNER1.out.progress_log )
    ch_log = ch_log.mix( CHECKM2_BINNER2.out.progress_log )
    ch_log = ch_log.mix( CHECKM2_BINNER3.out.progress_log )
    ch_log = ch_log.mix( CONSOLIDATE_BINS.out.progress_log )

    emit:
    bin_ref_bins = CONSOLIDATE_BINS.out.dereplicated_bins
    versions     = ch_versions
    progress_log = ch_log
}
