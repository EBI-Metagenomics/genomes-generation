/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP } from '../modules/utils'
include { CHANGE_DOT_TO_UNDERSCORE } from '../modules/utils'
include { CHANGE_DOT_TO_UNDERSCORE_READS } from '../modules/utils'
include { FASTP } from '../modules/fastp'
include { DECONTAMINATION } from '../subworkflows/decontamination'
include { CHANGE_ERR_TO_ERZ } from '../modules/utils'

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
    CHANGE_ERR_TO_ERZ(input_data.map{item -> tuple(item[0], item[2])}, rename_file.first())              // tuple(accession, [reads]])
    CHANGE_DOT_TO_UNDERSCORE_READS(CHANGE_ERR_TO_ERZ.out.return_files)

    // --- trimming reads
    FASTP(CHANGE_DOT_TO_UNDERSCORE_READS.out.underscore_reads)

    DECONTAMINATION(FASTP.out.output_reads, ref_genome.first(), ref_genome_name)
    emit:
        contigs_fixed = CHANGE_DOT_TO_UNDERSCORE.out.underscore_contigs
        reads_cleaned = DECONTAMINATION.out.decontaminated_reads
}
