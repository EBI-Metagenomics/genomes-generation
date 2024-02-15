process CONSOLIDATE_BINS {

    label 'process_low'
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
    tuple val(meta), path("consolidated_bins"),                      emit: consolidated_bins
    tuple val(meta), path("consolidated_stats.tsv"), optional: true, emit: consolidated_stats
    tuple val(meta), path("dereplicated_bins/*"),    optional: true, emit: dereplicated_bins
    tuple val(meta), path("dereplicated_list.tsv") , optional: true, emit: dereplicated_list
    path "versions.yml"                            ,                 emit: versions
    path "progress.log"                            ,                 emit: progress_log

    script:
    """
    consolidate_bins.py -i binner1 binner2 binner3 binner12 binner13 binner23 binner123 -s stats -v

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
        numpy: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('numpy').version)")
    END_VERSIONS

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}
        binner1: \$(ls binner1 | wc -l), binner2: \$(ls binner2 | wc -l), binner3: \$(ls binner3 | wc -l), binner12: \$(ls binner12 | wc -l), binner23: \$(ls binner23 | wc -l), binner13: \$(ls binner13 | wc -l), binner123: \$(ls binner123 | wc -l), consolidated: \$(ls consolidated_bins | wc -l)
    END_LOGGING
    """
}