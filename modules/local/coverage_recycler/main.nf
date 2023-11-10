process COVERAGE_RECYCLER {

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(genomes, stageAs: "genomes_dir/*")
    path metabat_depth

    output:
    tuple val(meta), path("coverage/*_coverage")        , emit: coverage_dir
    tuple val(meta), path("coverage/*_contigs2bins.txt"), emit: coverage_contigs2bins
    path "versions.yml"                                 , emit: versions

    script:
    """
    zcat ${metabat_depth} > depth.txt
    cov_recycler.py -g genomes_dir -m depth.txt -n ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}