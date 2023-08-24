process BREADTH_DEPTH {

    tag "${meta.id} ${mag}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cmseq:1.0.4--pyhb7b1952_0':
        'quay.io/biocontainers/cmseq:1.0.4--pyhb7b1952_0' }"

    input:
    tuple val(meta), path(mag), path(bam), path(bai)

    output:
    tuple val(meta), path("${mag.baseName}.coverage.csv"), emit: coverage

    script:
    """
    breadth_depth.py \
          --combine \
          --mincov 1 \
          ${bam} > "${mag.baseName}.coverage.csv"
    """

    stub:
    """
    touch \
        "${mag.baseName}.coverage.csv"
    """
}