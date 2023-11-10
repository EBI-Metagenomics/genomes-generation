process BREADTH_DEPTH {

    tag "${meta.id} ${mag}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cmseq:1.0.4--pyhb7b1952_0':
        'quay.io/biocontainers/cmseq:1.0.4--pyhb7b1952_0' }"

    input:
    tuple val(meta), path(mag), path(bam), path(bai)

    output:
    tuple val(meta), path("${mag.baseName}.coverage.csv"), emit: coverage
    path "versions.yml"                                  , emit: versions

    script:
    """
    breadth_depth.py \
          --combine \
          --mincov 1 \
          ${bam} > "${mag.baseName}.coverage.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """

    stub:
    """
    touch \
        "${mag.baseName}.coverage.csv"
    """
}