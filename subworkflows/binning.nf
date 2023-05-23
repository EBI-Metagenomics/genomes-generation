/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP } from '../modules/gunzip'
include { METAWRAP_BINNING } from '../modules/metawrap'
include { BIN_REFINEMENT } from '../modules/metawrap'

workflow BINNING {
    take:
        contigs  // tuple(name, contigs)
        reads    // tuple(name, contigs)
    main:

    GUNZIP(reads.map{item -> item[1]})

    METAWRAP_BINNING(contigs, GUNZIP.out.uncompressed)

    BIN_REFINEMENT(
        METAWRAP_BINNING.out.binning_metabat2,
        METAWRAP_BINNING.out.binning_concoct,
        METAWRAP_BINNING.out.binning_maxbin2)

    emit:
        concoct_bins = METAWRAP_BINNING.out.binning_concoct
        metabat2_bins = METAWRAP_BINNING.out.binning_metabat2
        binning_result = BIN_REFINEMENT.out.bin_ref_bins
        metabat_depth_for_coverage = METAWRAP_BINNING.out.metabat_depth_for_coverage
}
