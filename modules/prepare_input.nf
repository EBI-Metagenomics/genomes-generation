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
    path single_reads
    path paired_forward
    path paired_reverse

    output:
    path 'cleaning_reads/${name}_trimmed.fq*', emit: single_trimmed, optional: true
    path 'cleaning_reads/${name}_1_val_1.fq*', emit: paired_forward_trimmed, optional: true
    path 'cleaning_reads/${name}_2_val_2.fq*', emit: paired_reverse_trimmed, optional: true

    script:
    if (mode == 'paired') {
        """
        echo '[ fastq clean-up ] Cleaning FASTQ files'
        trim_galore --paired ${paired_forward} ${paired_reverse} -o cleaning_reads
        """
    }
    else if (mode == 'single') {
        """
        echo '[ fastq clean-up ] Cleaning FASTQ files'
        trim_galore ${single_reads} -o cleaning_reads
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
    path single_reads
    path paired_forward
    path paired_reverse

    output:
    path '${name}_clean.fastq.gz', emit: single_cleaned, optional: true
    path '${name}_clean_1.fastq.gz', emit: paired_forward_cleaned, optional: true
    path '${name}_clean_2.fastq.gz', emit: paired_reverse_cleaned, optional: true

    script:
    if (mode == 'paired') {
        """
        echo '[ fastq clean-up ] Mapping files to host genome'
        bwa-mem2 mem -M -t ${task.cpus} ${ref_genome}/${ref_genome_name} ${paired_forward} ${paired_reverse} | samtools view -@ ${task.cpus} -f 12 -F 256 -uS - -o ${name}_both_unmapped.bam
        samtools sort -@ ${task.cpus} -n ${name}_both_unmapped.bam -o ${name}_both_unmapped_sorted.bam
	    bedtools bamtofastq -i ${name}_both_unmapped_sorted.bam -fq ${name}_clean_1.fastq -fq2 ${name}_clean_2.fastq

	    echo '[ fastq clean-up ] Compressing output files'
	    gzip ${name}_clean_1.fastq
	    gzip ${name}_clean_2.fastq
        """
    }
    else if (mode == 'single') {
        """
        echo '[ fastq clean-up ] Mapping files to host genome'
        bwa-mem2 mem -M -t ${task.cpus} ${ref_genome}/${ref_genome_name} ${single_reads} | samtools view -@ ${task.cpus} -f 4 -F 256 -uS - -o ${name}_unmapped.bam
        samtools sort -@ ${task.cpus} -n ${name}_unmapped.bam -o ${name}_unmapped_sorted.bam
        bedtools bamtofastq -i ${name}_unmapped_sorted.bam -fq ${name}_clean.fastq

	    echo '[ fastq clean-up ] Compressing output file'
	    gzip ${name}_clean.fastq
        """
    }
}