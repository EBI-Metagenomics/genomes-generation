process BINNING_REFINER {

    label 'process_low'
    tag "${name} ${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    val(name)
    tuple val(meta), path(bin1, stageAs: "binner1/*"), path(bin2, stageAs: "binner2/*"), path(bin3, stageAs: "binner3/*")

    output:
    tuple val(meta), path("${meta.id}_output_${name}/refined"), emit: refined_bins
    path "versions.yml"                                       , emit: versions
    path "progress.log"                                       , emit: progress_log

    script:
    """
    binning_refiner.py -1 binner1 -2 binner2 -3 binner3 -o "${meta.id}_output_${name}" -n "${meta.id}_${name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}\t${name}
        binner1: \$(ls binner1 | wc -l), binner2: \$(ls binner2 | wc -l), binner3: \$(ls binner3 | wc -l), refined: \$(ls ${meta.id}_output_${name}/refined | wc -l)
    END_LOGGING
    """
}