/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT } from '../modules/align_bwa'
include { BAM_TO_FASTQ } from '../modules/bam_to_fastq'

workflow DECONTAMINATION {
    take:
        reads    // tuple(name, contigs)
        ref_genome
        ref_genome_name

    main:
    // TODO: fix mode
    def mode = "paired"
    //if (reads_list.size() == 1) {
    //    mode = 'single' }
    //else if (reads_list.size() == 2) {
    //    mode = 'paired' }
    //else {
    //    print('incorrect reads') }

    ALIGNMENT(reads, ref_genome, ref_genome_name, channel.value("-f 12 -F 256 -uS"))

    BAM_TO_FASTQ(ALIGNMENT.out.bams, mode)

    emit:
        decontaminated_reads = BAM_TO_FASTQ.out.clean_reads
}