process COVERAGE_RECYCLER {

    tag "${meta.id} ${metabat_depth}"

    publishDir(
        path: "${params.outdir}/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    tuple val(meta), path(genomes, stageAs: "genomes_dir/*")
    path(metabat_depth)

    output:
    tuple val(meta), path("coverage/*_coverage"), emit: coverage_dir
    tuple val(meta), path("coverage/*_contigs2bins.txt"), emit: coverage_contigs2bins

    script:
    """
    cov_recycler.py -g genomes_dir -m ${metabat_depth} -n ${meta.id}
    """
}