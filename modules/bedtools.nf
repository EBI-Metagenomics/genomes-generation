process BEDTOOLS_BAMTOFASTQ {

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
    val mode
    val name
    path bam

    output:
    path '*_clean*.fastq.gz', emit: reads_cleaned

    script:
    if (mode == 'paired') {
        """
	    bedtools bamtofastq -i ${bam} -fq ${name}_clean_1.fastq -fq2 ${name}_clean_2.fastq
	    gzip ${name}_clean_1.fastq
	    gzip ${name}_clean_2.fastq
        """
    }
    else if (mode == 'single') {
        """
        bedtools bamtofastq -i ${bam} -fq ${name}_clean.fastq
        gzip ${name}_clean.fastq
        """
    }
}