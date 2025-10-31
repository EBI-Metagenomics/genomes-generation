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
    tuple val(meta), path("*{_1,_2,_interleaved}.decont.fq.gz"), emit: reads
    path "versions.yml",              emit: versions

    script:
    def bam2fq_args = task.ext.bam2fq_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_args = task.ext.decontamination_args

    """
    if [[ "${meta.single_end}" == "true" ]]; then
        bwa-mem2 \
            mem \
            -M \
            -t $task.cpus \
            $ref_fasta \
            $reads \
            | samtools view -@ ${task.cpus} ${samtools_args} - \
            | samtools sort -@ ${task.cpus} -n -O bam - \
            | samtools bam2fq ${bam2fq_args} -@ $task.cpus - | gzip --no-name > ${prefix}_interleaved.decont.fq.gz
    else
        bwa-mem2 \
            mem \
            -M \
            -t $task.cpus \
            $ref_fasta \
            $reads \
            | samtools view -@ ${task.cpus} ${samtools_args} - \
            | samtools sort -@ ${task.cpus} -n -O bam - \
            | samtools bam2fq ${bam2fq_args} -@ ${task.cpus} -1 ${prefix}_1.decont.fq.gz -2 ${prefix}_2.decont.fq.gz -0 /dev/null -s /dev/null
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}