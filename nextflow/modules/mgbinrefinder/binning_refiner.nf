process BINNING_REFINER {

    tag "${name} ${meta.id}"

    container 'quay.io/biocontainers/biopython:1.75'

    publishDir(
        path: "${params.outdir}/intermediate_steps/Refinement/binref_${name}",
        mode: params.publish_dir_mode,
        failOnError: true,
        pattern: "output_${name}/Refined/*.fa"
    )

    input:
    val(name)
    tuple val(meta), path(bin1, stageAs: "binner1/*")
    tuple val(meta), path(bin2, stageAs: "binner2/*")

    output:
    tuple val(meta), path("output_${name}/Refined/*"), emit: refined_bins

    script:
    """
    binning_refiner.py -1 binner1 -2 binner2 -o output_${name} -n ${name}
    """
}

process BINNING_REFINER3 {

    tag "${name} ${meta.id}"

    container 'quay.io/biocontainers/biopython:1.75'

    publishDir(
        path: "${params.outdir}/binref_${name}",
        mode: 'copy',
        failOnError: true
    )

    input:
    val(name)
    tuple val(meta), path(bin1, stageAs: "binner1/*")
    tuple val(meta), path(bin2, stageAs: "binner2/*")
    tuple val(meta), path(bin3, stageAs: "binner3/*")

    output:
    tuple val(meta), path("output_${name}/Refined/*"), emit: refined_bins

    script:
    """
    binning_refiner.py -1 binner1 -2 binner2 -3 binner3 -o output_${name} -n ${name}
    """
}