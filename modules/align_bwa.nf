/*
 * Host decontamination
*/
process ALIGNMENT {

    label 'alignment'

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(name), path(input_ch_reads)
    path ref_genome
    val ref_genome_name

    output:
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
        ${reads} | samtools view -@ ${task.cpus} -f 4 -F 256 -uS - -o output_decontamination/${name}_unmapped.bam

        echo "samtools sort"
        samtools sort -@ ${task.cpus} -n output_decontamination/${name}_unmapped.bam \
        -o output_decontamination/${name}_unmapped_sorted.bam

        echo "samtools index sorted bam"
        samtools index output_decontamination/${name}_unmapped_sorted.bam
        """
    } else if ( reads.size() == 2 ) {
        if (reads[0].name.contains("_1")) {
            input_reads = "${reads[0]} ${reads[1]}"
        } else {
            input_reads = "${reads[1]} ${reads[0]}"
        }
        """
        mkdir -p output_decontamination
        echo "mapping files to host genome PE"
        bwa-mem2 mem -M \
        -t ${task.cpus} \
        ${ref_genome}/${ref_genome_name} \
        ${input_reads} | samtools view -@ ${task.cpus} -f 12 -F 256 -uS - -o output_decontamination/${name}_both_unmapped.bam

        echo "samtools sort"
        samtools sort -@ ${task.cpus} -n output_decontamination/${name}_both_unmapped.bam -o output_decontamination/${name}_both_unmapped_sorted.bam

        echo "samtools index sorted bam"
        samtools index output_decontamination/${name}_both_unmapped_sorted.bam
        """
    } else {
        error "Invalid mode input"
    }
}


process INDEX_REF_GENOME {

    label 'alignment'

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(name), path(genome)

    output:
    tuple val(name), path("ref_genome") , emit: index_folder

    script:
    """
    echo "Create index for ref genome"
    mkdir ref_genome
    cp ${genome} ref_genome
    cd ref_genome
    bwa-mem2 index ${genome}
    """
}