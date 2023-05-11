process PARSE_CHECKM {

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    path bins_taxonomy
    path bins_qa

    output:
    path "checkm_results.tab", emit: checkm_results

    script:
    """
    parse_checkm.py ${bins_qa} ${bins_taxonomy} > checkm_results.tab
    """
}