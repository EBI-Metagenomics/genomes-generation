process COVERAGE_RECYCLER {

    tag "${name}"

    publishDir(
        path: "${params.outdir}/coverage/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    tuple val(name), path(genomes, stageAs: "genomes_dir/*"), path(metabat_depth)

    output:
    tuple val(name), path("coverage"), emit: coverage_dir

    script:
    """
    cov_recycler.py -g genomes_dir -m ${metabat_depth}
    """
}