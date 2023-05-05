process SELECT_SEQS {

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    val sample_name
    path bins
    path cont_contigs

    output:
    path "_clean.fa", emit: clean_bins

    script:
    """
    python3 select_seqs_notin_IDs.py -i ${bins} -d ${cont_contigs} -o "${sample_name}_clean.fa"
    """
}