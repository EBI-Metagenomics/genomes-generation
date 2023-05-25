/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP } from '../modules/gunzip'
include { CHANGE_DOT_TO_UNDERSCORE } from '../modules/prepare_data'
include { FASTP } from '../modules/fastp'
include { DECONTAMINATION } from '../modules/decontamination'
include { CHANGE_ERR_TO_ERZ } from '../modules/prepare_data'

workflow PREPARE_INPUT {
    take:
        input_data
        ref_genome
        ref_genome_name
        rename_file
    main:
    // --- fix contig names
    // GUNZIP(name, contigs)
    CHANGE_DOT_TO_UNDERSCORE(input_data.map{item -> tuple(item[0], item[1])})           // tuple(accession, assembly)

    // change ERR in reads to ERZ
    CHANGE_ERR_TO_ERZ(input_data.map{item -> tuple(item[0], item[2])}, rename_file)              // tuple(accession, [reads]])

    // --- trimming reads
    FASTP(CHANGE_ERR_TO_ERZ.out.return_files)

    DECONTAMINATION(FASTP.out.output_reads, ref_genome.first(), ref_genome_name)
    emit:
        contigs_fixed = CHANGE_DOT_TO_UNDERSCORE.out.return_contigs
        reads_cleaned = DECONTAMINATION.out.decontaminated_reads
        bams = DECONTAMINATION.out.bams
}
