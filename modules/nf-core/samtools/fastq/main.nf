// FIXME: remove this module and install -> bam2fq
process SAMTOOLS_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(inputbam)
    val split

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (split){
        """
        samtools \\
            bam2fq \\
            $args \\
            -@ $task.cpus \\
            -1 ${prefix}_1.fastq.gz \\
            -2 ${prefix}_2.fastq.gz \\
            -0 ${prefix}_other.fastq.gz \\
            -s ${prefix}_singleton.fastq.gz \\
            $inputbam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        samtools \\
            bam2fq \\
            $args \\
            -@ $task.cpus \\
            $inputbam | gzip --no-name > ${prefix}_interleaved.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
