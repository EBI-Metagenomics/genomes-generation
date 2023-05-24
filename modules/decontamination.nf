/*
 * Host decontamination
*/
process DECONTAMINATION {

    publishDir "${params.outdir}/qc/decontamination", mode: 'copy'

    label 'decontamination'

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(name), path(input_ch_reads)
    path ref_genome
    val ref_genome_name

    output:
    tuple val(name), path("*_clean*.fastq.gz") , emit: decontaminated_reads
    tuple val(name), path("output_decontamination/${name}*_unmapped_sorted.bam*") , emit: bams

    script:
    reads = input_ch_reads.collect()
    def input_reads = "";

    if (reads.size() == 1 ) {
        input_reads = "${reads[0]}";
        """
        mkdir -p output_decontamination

        echo "mapping files to host genome SE"
        bwa-mem2 mem -M -t ${task.cpus} \
        ${ref_genome}/${ref_genome_name} \
        ${reads} > out.sam

        echo "convert sam to bam"
        samtools view -@ ${task.cpus} -f 4 -F 256 -uS -o output_decontamination/${name}_unmapped.bam out.sam

        echo "samtools sort"
        samtools sort -@ ${task.cpus} -n output_decontamination/${name}_unmapped.bam \
        -o output_decontamination/${name}_unmapped_sorted.bam

        echo "samtools index sorted bam"
        samtools index output_decontamination/${name}_unmapped_sorted.bam

        echo "samtools"
        samtools fastq output_decontamination/${name}_unmapped_sorted.bam > output_decontamination/${name}_clean.fastq

        echo "compressing output file"
        gzip -c output_decontamination/${name}_clean.fastq > ${name}_clean.fastq.gz
        """
    } else if ( reads.size() == 2 ) {
        if (reads[0].name.contains("_1")) {
            input_reads = "${reads[0]} ${reads[1]}"
        } else {
            input_reads = "${reads[1]} ${reads[0]}"
        }
        """
        mkdir output_decontamination
        echo "mapping files to host genome PE"
        bwa-mem2 mem -M \
        -t ${task.cpus} \
        ${ref_genome}/${ref_genome_name} \
        ${input_reads} > out.sam

        echo "convert sam to bam"
        samtools view -@ ${task.cpus} -f 12 -F 256 -uS -o output_decontamination/${name}_both_unmapped.bam out.sam

        echo "samtools sort"
        samtools sort -@ ${task.cpus} -n output_decontamination/${name}_both_unmapped.bam -o output_decontamination/${name}_both_unmapped_sorted.bam

        echo "samtools index sorted bam"
        samtools index output_decontamination/${name}_both_unmapped_sorted.bam

        echo "samtools fastq"
        samtools fastq -1 output_decontamination/${name}_clean_1.fastq \
        -2 output_decontamination/${name}_clean_2.fastq \
        -0 /dev/null \
        -s /dev/null \
        -n output_decontamination/${name}_both_unmapped_sorted.bam

        echo "compressing output files"
        gzip -c output_decontamination/${name}_clean_1.fastq > ${name}_clean_1.fastq.gz
        gzip -c output_decontamination/${name}_clean_2.fastq > ${name}_clean_2.fastq.gz
        """
    } else {
        error "Invalid mode input"
    }
}