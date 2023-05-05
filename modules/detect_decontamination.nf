process DETECT_DECONTAMINATION {

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    val sample_name
    path summary
    path names

    output:
    path "*.cont-contigs.txt", emit: cont_contigs
    path "*.cont-stats.tsv", emit: cont_stats

    script:
    """
    python3 detect_contamination.py -s ${summary} -n ${names} -i ${sample_name}
    """
}