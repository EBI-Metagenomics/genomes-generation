/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_DOT_TO_UNDERSCORE {

    container 'quay.io/biocontainers/drep:3.2.2--pyhdfd78af_0'

    input:
    path contigs

    output:
    path "out", emit: contigs

    script:
    """
    bash sed.sh ${contigs}
    """
}

/*
 * trimming
*/
process TRIM_GALORE {

    container 'quay.io/microbiome-informatics/trim_galore:0.6.7'

    input:
    val mode
    val name
    path reads

    output:
    path 'cleaning_reads/*.fq', emit: reads_trimmed

    script:
    if (mode == 'paired') {
        """
        echo '[ fastq clean-up ] Cleaning FASTQ files'
        trim_galore --paired ${reads[0]} ${reads[1]} -o cleaning_reads
        """
    }
    else if (mode == 'single') {
        """
        echo '[ fastq clean-up ] Cleaning FASTQ files'
        trim_galore ${reads} -o cleaning_reads
        """
    }
}

process MAP_HOST_GENOME {

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    val mode
    val name
    path ref_genome
    val ref_genome_name
    path reads

    output:
    path '*_unmapped_sorted.bam', emit: bam_sorted

    script:
    if (mode == 'paired') {
        """
        echo '[ fastq clean-up ] Mapping files to host genome'
        bwa-mem2 mem -M -t ${task.cpus} ${ref_genome}/${ref_genome_name} ${reads[0]} ${reads[1]} | samtools view -@ ${task.cpus} -f 12 -F 256 -uS - -o ${name}_both_unmapped.bam
        echo '[ fastq clean-up ] sort bam'
        samtools sort -@ ${task.cpus} -n ${name}_both_unmapped.bam -o ${name}_both_unmapped_sorted.bam
        """
    }
    else if (mode == 'single') {
        """
        echo '[ fastq clean-up ] Mapping files to host genome'
        bwa-mem2 mem -M -t ${task.cpus} ${ref_genome}/${ref_genome_name} ${reads} | samtools view -@ ${task.cpus} -f 4 -F 256 -uS - -o ${name}_unmapped.bam
        echo '[ fastq clean-up ] sort bam'
        samtools sort -@ ${task.cpus} -n ${name}_unmapped.bam -o ${name}_unmapped_sorted.bam
        """
    }
}