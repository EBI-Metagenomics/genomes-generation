process DETECT_DECONTAMINATION {

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    tuple val(sample_name), path(summary)
    tuple val(sample_name), path(names)

    output:
    tuple val(sample_name), path("*.cont-contigs.txt"), emit: cont_contigs
    tuple val(sample_name), path("*.cont-stats.tsv"), emit: cont_stats

    script:
    """
    detect_decontamination.py -s ${summary} -n ${names} -i ${sample_name}
    """
}