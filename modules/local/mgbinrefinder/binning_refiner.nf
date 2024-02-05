process BINNING_REFINER {

    label 'process_low'
    tag "${name} ${meta.id}"

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    val(name)
    tuple val(meta), path(bin1, stageAs: "binner1/*")
    tuple val(meta), path(bin2, stageAs: "binner2/*")
    tuple val(meta), path(bin3, stageAs: "binner3/*")

    output:
    tuple val(meta), path("${meta.id}_output_${name}/refined/*"), optional: true, emit: refined_bins
    path "versions.yml"                                                         , emit: versions

    script:
    """
    binning_refiner.py -1 binner1 -2 binner2 -3 binner3 -o "${meta.id}_output_${name}" -n "${meta.id}_${name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}