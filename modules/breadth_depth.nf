process BREADTH_DEPTH {

    tag "${name}"
    container 'quay.io/biocontainers/cmseq:1.0.4--pyhb7b1952_0'

    publishDir(
        path: "${params.outdir}/euk_coverage/",
        mode: 'copy',
        failOnError: true
    )

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}.coverage.csv"), emit: coverage

    script:
    bam_files = bam.collect()
    def bam_file = ""
    if (bam_files[0].contains('.bai')) {
        bam_file = bam_files[1]
    }
    else {
        bam_file = bam_files[0]
    }
    """
    breadth_depth.py \
          --combine \
          --mincov 1 \
          ${bam_file} > "${name}.coverage.csv"
    """

    stub:
    """
    touch \
        "${name}.coverage.csv"
    """
}