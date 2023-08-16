process BREADTH_DEPTH {

    tag "${meta.id} ${mag}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cmseq:1.0.4--pyhb7b1952_0':
        'quay.io/biocontainers/cmseq:1.0.4--pyhb7b1952_0' }"

    input:
    tuple val(meta), path(bam), path(mag)

    output:
    tuple val(meta), path("${mag.baseName}.coverage.csv"), emit: coverage

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
          ${bam_file} > "${mag.baseName}.coverage.csv"
    """

    stub:
    """
    touch \
        "${mag.baseName}.coverage.csv"
    """
}