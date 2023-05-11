process COVERAGE_RECYCLER {

    publishDir(
        path: "${params.outdir}/coverage/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    val name
    path dereplicated_genomes
    path metabat_depth

    output:
    path "coverage", emit: coverage_dir

    script:
    """
    cov_recycler.py -g ${dereplicated_genomes} -m ${metabat_depth}
    """
}