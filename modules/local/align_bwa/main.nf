
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

process ALIGNMENT_DEPTH {

    /*
    This module aligns reads to the reference using specified arguments and produces output in the form of FASTQ.GZ, BAM and BAM.BAI files.
    metabat2 jgi_summarize_bam_contig_depths generates depths of input reads.
    BAM files removed in the end.
    */

    label 'process_medium'

    tag "${meta.id} align to ${ref_fasta}"

    container 'quay.io/microbiome-informatics/bwa_metabat:2.2.1_2.15'

    input:
    tuple val(meta), path(reads), path(ref_fasta), path(ref_fasta_index)

    output:
    tuple val(meta), path(ref_fasta), path("*.txt.gz"), emit: depth
    path "versions.yml"                               , emit: versions

    script:
    def input_reads = "";
    def args_jgi_summarize_bam_contig_depths = task.ext.args_jgi_summarize_bam_contig_depths ?: ''
    def prefix_jgi_summarize_bam_contig_depths = task.ext.prefix_jgi_summarize_bam_contig_depths ?: "${meta.id}"

    if ( meta.single_end ) {
        input_reads = "${reads[0]}";
    } else {
        input_reads = "${reads[0]} ${reads[1]}"
    }

    def samtools_args = task.ext.alignment_args

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

    echo "generate depth"
    export OMP_NUM_THREADS=$task.cpus

    jgi_summarize_bam_contig_depths \\
        --outputDepth ${prefix_jgi_summarize_bam_contig_depths}.txt \\
        $args_jgi_summarize_bam_contig_depths \\
        output/${meta.id}_sorted.bam

    echo "compress depth"
    bgzip --threads $task.cpus ${prefix_jgi_summarize_bam_contig_depths}.txt

    echo "remove bams"
    rm -rf output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS
    """
}

process ALIGNMENT_DECONTAMINATION {

    /*
    This module aligns reads to the reference using specified arguments and produces output in the form of FASTQ.GZ files.
    Once the alignment is completed, the BAM files are removed to save disk space.
    */

    label 'process_medium'

    tag "${meta.id} align to ${ref_fasta}"

    container 'quay.io/microbiome-informatics/bwamem2:2.2.1'

    input:
    tuple val(meta), path(reads), path(ref_fasta), path(ref_fasta_index)

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

    def samtools_args = task.ext.decontamination_args

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

    echo "remove bams"
    rm -rf output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}