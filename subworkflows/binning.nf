/*
    ~~~~~~~~~~~~~~~~~~
     Collect bins
    ~~~~~~~~~~~~~~~~~~
*/
process COLLECT_BINS {
    container 'quay.io/microbiome-informatics/metawrap:latest'

    input:
    val name
    path bin_ref_folder

    output:
    path "Metawrap_bins", emit: metawrap_bins

    script:
    """
    mkdir Metawrap_bins

    """
}

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
        mode
        name
        contigs
        reads
    main:

    reads_list = reads.collect()
    GUNZIP(name, reads_list)

    METAWRAP_BINNING(mode, name, contigs, GUNZIP.out.uncompressed)

    BIN_REFINEMENT(
        METAWRAP_BINNING.out.binning_metabat2,
        METAWRAP_BINNING.out.binning_concoct,
        METAWRAP_BINNING.out.binning_maxbin2)

    //COLLECT_BINS(name, BIN_REFINEMENT.out.bin_ref)

    emit:
        binning_result = BIN_REFINEMENT.out.bin_ref
}
