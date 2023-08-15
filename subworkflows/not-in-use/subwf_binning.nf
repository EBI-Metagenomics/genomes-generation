/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP } from '../modules/utils'
include { METAWRAP_BINNING as METAWRAP_BINNING_METABAT2} from '../modules/metawrap'
include { METAWRAP_BINNING as METAWRAP_BINNING_CONCOCT} from '../modules/metawrap'
include { METAWRAP_BINNING as METAWRAP_BINNING_MAXBIN2} from '../modules/metawrap'

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

    // construct output
    // tuple( run_accession, assembly_file, [raw_reads], concoct_folder, metabat_folder )
    // euk
    output_for_euk_part = input_data.combine(METAWRAP_BINNING_CONCOCT.out.binning, by: 0).combine(METAWRAP_BINNING_METABAT2.out.binning, by: 0)
    // prok
    collect_binners = METAWRAP_BINNING_METABAT2.out.binning.combine(METAWRAP_BINNING_CONCOCT.out.binning, by: 0).combine(METAWRAP_BINNING_MAXBIN2.out.binning, by: 0)
    metabat_depth = METAWRAP_BINNING_METABAT2.out.metabat_depth_for_coverage.collectFile(name:"aggregated_metabat_depth.txt", skip:1, keepHeader:true)

    emit:
        collect_binners
        output_for_euk_part
        metabat_depth
}
