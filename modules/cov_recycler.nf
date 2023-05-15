process COVERAGE_RECYCLER {

    publishDir(
        path: "${params.outdir}/coverage/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    path genomes, stageAs: "genomes_dir/*"
    path metabat_depth

    output:
    path "coverage", emit: coverage_dir

    script:
    """
    cov_recycler.py -g genomes_dir -m ${metabat_depth}
    """
}