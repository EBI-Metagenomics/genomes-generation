include { BINNING_REFINER   } from '../../modules/mgbinrefinder/binning_refiner'
include { CHECKM2           } from '../../modules/checkm2'

workflow REFINE {
    take:
        name
        binner1
        binner2
        binner3
        checkm_db
    main:
        BINNING_REFINER(name, binner1, binner2, binner3)
        size_refined_bins = BINNING_REFINER.out.refined_bins.collect().size()
        size_refined_bins.subscribe { println "Refinder: $name.value: $it.value" }

        CHECKM2(name, BINNING_REFINER.out.refined_bins, checkm_db)
        //size_filtered_bins = CHECKM2.out.checkm2_results_filtered.collect().size()
        //size_filtered_bins.subscribe { println "Checkm2: $name.value: $it.value" }
    emit:
        refined = BINNING_REFINER.out.refined_bins
        filtered_bins = CHECKM2.out.checkm2_results_filtered
        filtered_bins_stats = CHECKM2.out.checkm2_results_filtered_stats
}
