/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP } from '../modules/utils'
include { METAWRAP_BINNING as METAWRAP_BINNING_METABAT2} from '../modules/metawrap'
include { METAWRAP_BINNING as METAWRAP_BINNING_CONCOCT} from '../modules/metawrap'
include { METAWRAP_BINNING as METAWRAP_BINNING_MAXBIN2} from '../modules/metawrap'
include { BIN_REFINEMENT } from '../modules/metawrap'

workflow BINNING {
    take:
        accession_reads             // tuple(accession, [reads])
        accession_contigs           // tuple(accession, contigs)
    main:

    GUNZIP(accession_reads.map{item -> item[1]})            // reads

    METAWRAP_BINNING_METABAT2(accession_contigs, GUNZIP.out.uncompressed, channel.value("metabat2"))  // tuple(name, contigs), reads
    METAWRAP_BINNING_CONCOCT(accession_contigs, GUNZIP.out.uncompressed, channel.value("concoct"))
    METAWRAP_BINNING_MAXBIN2(accession_contigs, GUNZIP.out.uncompressed, channel.value("maxbin2"))

    BIN_REFINEMENT(
        METAWRAP_BINNING_METABAT2.out.binning,
        METAWRAP_BINNING_CONCOCT.out.binning,
        METAWRAP_BINNING_MAXBIN2.out.binning)

    emit:
        concoct_bins = METAWRAP_BINNING_CONCOCT.out.binning  // (name, folder)
        metabat2_bins = METAWRAP_BINNING_METABAT2.out.binning
        binning_result = BIN_REFINEMENT.out.bin_ref_bins
        metabat_depth_for_coverage = METAWRAP_BINNING_METABAT2.out.metabat_depth_for_coverage
}
