/*
 * Host decontamination
*/
process BAM_TO_FASTQ {
    tag "${name}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    label 'bam_to_fastq'

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(name), path(input_ch_bam)
    val(mode)

    output:
    tuple val(name), path("*_clean*.fastq.gz") , emit: clean_reads

    script:
    bams = input_ch_bam.collect()
    if ( mode == "single" ) {
        """
        echo "samtools"
        samtools fastq -@ ${task.cpus} ${bams[0]} > ${name}_clean.fastq

        gzip ${name}_clean.fastq
        """
    } else if ( mode == "paired" ) {
        """
        echo "samtools fastq"
        samtools fastq -@ ${task.cpus} \
        -1 ${name}_clean_1.fastq \
        -2 ${name}_clean_2.fastq \
        -0 /dev/null \
        -s /dev/null \
        -n ${bams[0]}

        gzip ${name}_clean_1.fastq ${name}_clean_2.fastq
        """
    } else {
        error "Invalid mode input"
    }
}