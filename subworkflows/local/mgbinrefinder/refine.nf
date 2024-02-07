include { BINNING_REFINER                } from '../../../modules/local/mgbinrefinder/binning_refiner'
include { CHECKM2 as CHECKM2_REFINE      } from '../../../modules/local/checkm2/main'

workflow REFINE {
    take:
    name
    binners
    checkm2_db

    main:

    ch_versions = Channel.empty()

    BINNING_REFINER( name, binners )
    ch_versions = ch_versions.mix( BINNING_REFINER.out.versions.first() )

    CHECKM2_REFINE( name, BINNING_REFINER.out.refined_bins, checkm2_db )
    ch_versions = ch_versions.mix( CHECKM2_REFINE.out.versions.first() )

    emit:
    refined = BINNING_REFINER.out.refined_bins
    filtered_bins = CHECKM2_REFINE.out.filtered_genomes
    filtered_bins_stats = CHECKM2_REFINE.out.filtered_stats
    versions = ch_versions
}


