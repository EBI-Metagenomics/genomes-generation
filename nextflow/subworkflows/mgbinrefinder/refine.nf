include { BINNING_REFINER                } from '../../modules/mgbinrefinder/binning_refiner'
include { BINNING_REFINER3               } from '../../modules/mgbinrefinder/binning_refiner'
include { CHECKM2 as CHECKM2_REFINE      } from '../../modules/local/checkm2/main'

workflow REFINE {
    take:
        name
        binner1
        binner2
        binner3
        checkm_db
    main:
        refined = Channel.empty()
        if (binner3) {
            BINNING_REFINER3(name, binner1, binner2, binner3)
            refined = BINNING_REFINER3.out.refined_bins
        }
        else {
            BINNING_REFINER(name, binner1, binner2)
            refined = BINNING_REFINER.out.refined_bins
        }

        CHECKM2_REFINE(name, refined.map{it -> [it[0], it[1].listFiles().flatten()]}, checkm_db.first())

    emit:
        refined = refined
        filtered_bins = CHECKM2_REFINE.out.checkm2_results_filtered
        filtered_bins_stats = CHECKM2_REFINE.out.checkm2_results_filtered_stats
}


