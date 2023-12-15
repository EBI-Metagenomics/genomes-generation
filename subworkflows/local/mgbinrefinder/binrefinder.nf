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

    REFINE12( "binner12", renamed_binner1, renamed_binner2, false, ref_checkm )
    REFINE13( "binner13", renamed_binner1, renamed_binner3, false, ref_checkm )
    REFINE23( "binner23", renamed_binner2, renamed_binner3, false, ref_checkm )
    REFINE123( "binner123", renamed_binner1, renamed_binner2, renamed_binner3, ref_checkm )

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

    CONSOLIDATE_BINS( binners, stats )
    empty_result = CONSOLIDATE_BINS.out.consolidated_stats.map{ meta, csv -> return tuple(meta, []) }

    ch_versions = ch_versions.mix( CONSOLIDATE_BINS.out.versions.first() )

    emit:
    bin_ref_bins = CONSOLIDATE_BINS.out.dereplicated_bins.ifEmpty(empty_result)
    versions     = ch_versions
}
