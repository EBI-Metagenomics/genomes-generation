/*
 * Host decontamination
*/
process ALIGNMENT {
    tag "${name}"
    label 'alignment'

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(name), path(input_ch_reads)
    path ref_genome
    val ref_genome_name
    val samtools_args

    output:
    tuple val(name), path("output/${name}_sorted.bam*") , emit: bams

    script:
    // define reads
    reads = input_ch_reads.collect()
    def input_reads = "";
    if (reads.size() == 1 ) {
        input_reads = "${reads[0]}";
    } else if ( reads.size() == 2 ) {
        if (reads[0].name.contains("_1")) {
            input_reads = "${reads[0]} ${reads[1]}"
        } else {
            input_reads = "${reads[1]} ${reads[0]}"
        }
    } else {
        error "Invalid mode input"
    }

    // align
    """
    mkdir -p output
    echo "mapping files to host genome"
    bwa-mem2 mem -M \
      -t ${task.cpus} \
      ${ref_genome}/${ref_genome_name} \
      ${input_reads} | \
    samtools view -@ ${task.cpus} ${samtools_args} - | \
    samtools sort -@ ${task.cpus} -O bam - -o output/${name}_sorted.bam

    echo "samtools index sorted bam"
    samtools index output/${name}_sorted.bam
    """
}

/*
 * Host decontamination
*/
process ALIGNMENT_WITH_INDEXING {

    tag "${name} vs ${ref_genome_name}"

    label 'alignment'

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(name), path(input_ch_reads)
    path ref_genome
    val ref_genome_name
    val samtools_args

    output:
    tuple val(name), path("output/${name}_sorted.bam*") , emit: bams

    script:

    // define reads
    reads = input_ch_reads.collect()
    def input_reads = "";
    if (reads.size() == 1 ) {
        input_reads = "${reads[0]}";
    } else if ( reads.size() == 2 ) {
        if (reads[0].name.contains("_1")) {
            input_reads = "${reads[0]} ${reads[1]}"
        } else {
            input_reads = "${reads[1]} ${reads[0]}"
        }
    } else {
        error "Invalid mode input"
    }

    // align
    """
    mkdir -p output
    echo "index ref genome"
    bwa-mem2 index ${ref_genome}

    echo "mapping files to host genome"
    bwa-mem2 mem -M \
      -t ${task.cpus} \
      ${ref_genome} \
      ${input_reads} | \
    samtools view -@ ${task.cpus} ${samtools_args} - | \
    samtools sort -@ ${task.cpus} -O bam - -o output/${name}_sorted.bam

    echo "samtools index sorted bam"
    samtools index -@ ${task.cpus} output/${name}_sorted.bam
    """
}

process INDEX_FASTA_META {

    tag "${meta.id} index ${fasta}"

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("${fasta.baseName}*.?*"), emit: fasta_with_index

    script:
    """
    echo "index ref genome"
    bwa-mem2 index ${fasta}
    """
}

process ALIGNMENT_META {

    tag "${meta.id} align to ${ref_fasta}"

    label 'alignment'

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(input_ch_reads)
    tuple val(meta), path(ref_fasta), path(ref_fasta_index)
    val samtools_args

    output:
    tuple val(meta), path(ref_fasta), path("output/${meta.id}_sorted.bam"), path("output/${meta.id}_sorted.bam.bai") , emit: bams

    script:

    // define reads
    reads = input_ch_reads.collect()
    def input_reads = "";
    if (reads.size() == 1 ) {
        input_reads = "${reads[0]}";
    } else if ( reads.size() == 2 ) {
        if (reads[0].name.contains("_1")) {
            input_reads = "${reads[0]} ${reads[1]}"
        } else {
            input_reads = "${reads[1]} ${reads[0]}"
        }
    } else {
        error "Invalid mode input"
    }

    // align
    """
    mkdir -p output

    echo "mapping files to host genome"
    bwa-mem2 mem -M \
      -t ${task.cpus} \
      ${ref_fasta} \
      ${input_reads} | \
    samtools view -@ ${task.cpus} ${samtools_args} - | \
    samtools sort -@ ${task.cpus} -O bam - -o output/${meta.id}_sorted.bam

    echo "samtools index sorted bam"
    samtools index -@ ${task.cpus} output/${meta.id}_sorted.bam
    """
}