include { BINNING_REFINER                } from '../../../modules/local/mgbinrefinder/binning_refiner'
include { BINNING_REFINER3               } from '../../../modules/local/mgbinrefinder/binning_refiner'
include { CHECKM2 as CHECKM2_REFINE      } from '../../../modules/local/checkm2/main'

workflow REFINE {
    take:
    name
    binner1
    binner2
    binner3
    checkm2_db

    main:
    refined = Channel.empty()
    if ( binner3 ) {
        BINNING_REFINER3( name, binner1.join( binner2 ).join( binner3 ) )
        refined = BINNING_REFINER3.out.refined_bins
    } else {
        BINNING_REFINER( name, binner1.join( binner2 ) )
        refined = BINNING_REFINER.out.refined_bins
    }

    CHECKM2_REFINE( name, refined, checkm2_db )

    emit:
    refined = refined
    filtered_bins = CHECKM2_REFINE.out.filtered_genomes
    filtered_bins_stats = CHECKM2_REFINE.out.filtered_stats
    // TODO: versions
}


