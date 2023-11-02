process COVERAGE_RECYCLER {

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(genomes, stageAs: "genomes_dir/*")
    path(metabat_depth)

    output:
    tuple val(meta), path("coverage/*_coverage"), emit: coverage_dir
    tuple val(meta), path("coverage/*_contigs2bins.txt"), emit: coverage_contigs2bins

    script:
    """
    zcat ${metabat_depth} > depth.txt
    cov_recycler.py -g genomes_dir -m depth.txt -n ${meta.id}
    """
}