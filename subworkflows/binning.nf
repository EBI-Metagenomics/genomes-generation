/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP } from '../modules/utils'
include { METAWRAP_BINNING } from '../modules/metawrap'
include { BIN_REFINEMENT } from '../modules/metawrap'

workflow BINNING {
    take:
        input_data  // tuple(name, reads, contigs)
    main:

    GUNZIP(input_data.map{item -> item[1]})  // reads

    METAWRAP_BINNING(input_data.map{item -> tuple(item[0], item[2])}, GUNZIP.out.uncompressed)  // tuple(name, contigs), reads

    BIN_REFINEMENT(
        METAWRAP_BINNING.out.binning_metabat2,
        METAWRAP_BINNING.out.binning_concoct,
        METAWRAP_BINNING.out.binning_maxbin2)

    emit:
        concoct_bins = METAWRAP_BINNING.out.binning_concoct  // (name, folder)
        metabat2_bins = METAWRAP_BINNING.out.binning_metabat2
        binning_result = BIN_REFINEMENT.out.bin_ref_bins
        metabat_depth_for_coverage = METAWRAP_BINNING.out.metabat_depth_for_coverage
}
