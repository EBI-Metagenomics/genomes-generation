process BINNING_REFINER {

    tag "${name}"

    container 'quay.io/biocontainers/biopython:1.75'

    publishDir(
        path: "${params.outdir}/binref_${name}",
        mode: 'copy',
        failOnError: true
    )

    input:
    val(name)
    path(bin1, stageAs: "binner1/*")
    path(bin2, stageAs: "binner2/*")
    path(bin3, stageAs: "binner3/*")

    output:
    path("output_${name}/Refined/*"), emit: refined_bins

    script:
    def args = ""
    if (!(bin3.toString().contains('NO_FILE'))) {
        args = "-3 binner3"
    }
    """
    binning_refiner.py -1 binner1 -2 binner2 ${args} -o output_${name} -n ${name}
    """
}