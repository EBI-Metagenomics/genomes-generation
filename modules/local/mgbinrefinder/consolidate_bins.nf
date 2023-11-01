process CONSOLIDATE_BINS {
    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), 
          path(bins1, stageAs: "binner1"),
          path(bins2, stageAs: "binner2"),
          path(bins3, stageAs: "binner3"),
          path(bins12, stageAs: "binner12"),
          path(bins13, stageAs: "binner13"),
          path(bins23, stageAs: "binner23"),
          path(bins123, stageAs: "binner123")
    tuple val(meta),
          path(stats1, stageAs: "stats/*"),
          path(stats2, stageAs: "stats/*"),
          path(stats3, stageAs: "stats/*"),
          path(stats12, stageAs: "stats/*"),
          path(stats13, stageAs: "stats/*"),
          path(stats23, stageAs: "stats/*"),
          path(stats123, stageAs: "stats/*")

    output:
    tuple val(meta), path("consolidated_bins"), emit: consolidated_bins
    tuple val(meta), path("consolidated_stats.tsv"), emit: consolidated_stats
    tuple val(meta), path("dereplicated_bins/*"), emit: dereplicated_bins

    script:
    """
    consolidate_bins.py -i binner1 binner2 binner3 binner12 binner13 binner23 binner123 -s stats -v
    """
}