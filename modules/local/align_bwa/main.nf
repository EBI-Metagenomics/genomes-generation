
process INDEX_FASTA {

    label 'process_medium'
    tag "${meta.id} index ${fasta}"

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("${fasta.baseName}*.?*"), emit: fasta_with_index
    path("versions.yml"),                                        emit: versions

    script:
    """
    echo "index ref genome"
    bwa-mem2 index ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
    END_VERSIONS
    """
}

process ALIGNMENT_BAM {
    // module does alignment with chosen arguments and generates BAM and BAM.BAI files as output
    label 'process_medium'

    tag "${meta.id} align to ${ref_fasta}"

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(reads), path(ref_fasta), path(ref_fasta_index)
    val align  // if true: align (include reads), else: decontaminate (exclude reads)

    output:
    tuple val(meta), path(ref_fasta), path("output/${meta.id}_sorted.bam"), path("output/${meta.id}_sorted.bam.bai"), emit: bam
    path "versions.yml"                                                                                             , emit: versions

    script:
    def input_reads = "";
    if ( meta.single_end ) {
        input_reads = "${reads[0]}";
    } else {
        input_reads = "${reads[0]} ${reads[1]}"
    }

    def samtools_args = ""
    if ( align ) {
        samtools_args = task.ext.alignment_args
    } else {
        samtools_args = task.ext.decontamination_args
    }
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process ALIGNMENT_READS {

    // module does alignment with chosen arguments and generates FASTQ.GZ files as output
    // that module was done to remove storage of heavy bam files

    label 'process_medium'

    tag "${meta.id} align to ${ref_fasta}"

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(reads), path(ref_fasta), path(ref_fasta_index)
    val align  // if true: align (include reads), else: decontaminate (exclude reads)

    output:
    tuple val(meta), path("*_*.fq.gz"), emit: reads
    path "versions.yml",              emit: versions

    script:
    def bam2fq_args = task.ext.bam2fq_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input_reads = "";
    if ( meta.single_end ) {
        input_reads = "${reads[0]}";
    } else {
        input_reads = "${reads[0]} ${reads[1]}"
    }

    def samtools_args = ""
    if ( align ) {
        samtools_args = task.ext.alignment_args
    } else {
        samtools_args = task.ext.decontamination_args
    }

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

    if [[ "${meta.single_end}" == "true" ]]; then
       samtools \
            bam2fq \
            $bam2fq_args \
            -@ $task.cpus \
            output/${meta.id}_sorted.bam | gzip --no-name > ${prefix}_interleaved.fq.gz
    else
        samtools \
            bam2fq \
            $bam2fq_args \
            -@ $task.cpus \
            -1 ${prefix}_1.fq.gz \
            -2 ${prefix}_2.fq.gz \
            -0 ${prefix}.other.fq.gz \
            -s ${prefix}.singleton.fq.gz \
            output/${meta.id}_sorted.bam
    fi

    rm -rf output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}