process CONSOLIDATE_BINS {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'biocontainers/biopython:1.75' }"

    input:
    path(binner_folders)
    path(stats_files, stageAs: "stats/*")

    output:
    path("consolidated_bins"), emit: consolidated_bins
    path("consolidated_stats.tsv"), emit: consolidated_stats
    path("dereplicated_bins/*"), emit: dereplicated_bins

    script:
    """
    consolidate_bins.py -i ${binner_folders} -s stats -v
    """
}