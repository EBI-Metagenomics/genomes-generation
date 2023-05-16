/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP } from '../modules/gunzip'
include { CHANGE_DOT_TO_UNDERSCORE } from '../modules/prepare_data'
include { FASTP } from '../modules/fastp'
include { DECONTAMINATION } from '../modules/decontamination'

workflow PREPARE_INPUT {
    take:
        input_data
        ref_genome
        ref_genome_name
    main:
    // --- fix contig names
    // GUNZIP(name, contigs)
    CHANGE_DOT_TO_UNDERSCORE(input_data.map{item -> tuple(item[0], item[1])})           // tuple(accession, assembly)

    // --- trimming reads
    FASTP(input_data.map{item -> tuple(item[0], item[2])})                              // tuple(accession, [reads]])
    DECONTAMINATION(FASTP.out.output_reads, ref_genome, ref_genome_name)

    emit:
        contigs_fixed = CHANGE_DOT_TO_UNDERSCORE.out.return_contigs
        reads_cleaned = DECONTAMINATION.out.decontaminated_reads
        bams = DECONTAMINATION.out.bams
}
