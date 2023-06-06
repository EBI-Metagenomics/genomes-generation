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
        input_data             // tuple( run_accession, assembly_file, [raw_reads] )
    main:

    reads = input_data.map(item -> item[2])
    accession_contigs = input_data.map(item -> tuple(item[0], item[1]))

    GUNZIP(reads)

    METAWRAP_BINNING_METABAT2(accession_contigs, GUNZIP.out.uncompressed, channel.value("metabat2"))  // tuple(name, contigs), reads
    METAWRAP_BINNING_CONCOCT(accession_contigs, GUNZIP.out.uncompressed, channel.value("concoct"))
    METAWRAP_BINNING_MAXBIN2(accession_contigs, GUNZIP.out.uncompressed, channel.value("maxbin2"))

    collect_binners = METAWRAP_BINNING_METABAT2.out.binning.combine(METAWRAP_BINNING_CONCOCT.out.binning, by: 0).combine(METAWRAP_BINNING_MAXBIN2.out.binning, by: 0)
    BIN_REFINEMENT(collect_binners)

    // construct output
    // tuple( run_accession, assembly_file, [raw_reads], concoct_folder, metabat_folder )
    output_for_euk_part = input_data.combine(METAWRAP_BINNING_CONCOCT.out.binning, by: 0).combine(METAWRAP_BINNING_METABAT2.out.binning, by: 0)

    // tuple( run_accession, bin_refinement, depth_file )
    output_for_prok_part = BIN_REFINEMENT.out.bin_ref_bins.combine(METAWRAP_BINNING_METABAT2.out.metabat_depth_for_coverage, by: 0)

    emit:
        output_for_euk_part
        output_for_prok_part
}
