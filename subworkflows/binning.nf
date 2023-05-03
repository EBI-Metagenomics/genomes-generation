include { METAWRAP_BINNING } from '../modules/metawrap'
include { BIN_REFINEMENT } from '../modules/metawrap'

process COLLECT_BINS {
    input:
        val name
        path bin_ref_folder
    output: "Metawrap_bins", emit metawrap_bins
    script:
    """
    mkdir Metawrap_bins
    for y in ${bin_ref_folder}
    do
        cp ${bin_ref_folder}/$y Metawrap_bins/${name}_$y
    done
    """
}

/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow BINNING {
    take:
        mode
        name
        contigs
        reads
    main:

    METAWRAP_BINNING(mode, name, contigs, reads)

    BIN_REFINEMENT(
        METAWRAP_BINNING.out.binning_metabat2,
        METAWRAP_BINNING.out.binning_concoct,
        METAWRAP_BINNING.out.binning_maxbin2)

    emit:
        binning_result = BIN_REFINEMENT.bin_ref
}
