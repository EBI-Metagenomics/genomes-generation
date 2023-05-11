process SELECT_SEQS {

    container 'quay.io/biocontainers/biopython:1.75'
    publishDir(
        path: "${params.outdir}/cleaned_bins",
        mode: "copy"
    )

    input:
    val sample_name
    path bins
    path cont_contigs

    output:
    path "*_clean.fa", emit: clean_bins

    script:
    """
    select_seqs_notin_IDs.py -i ${bins} -d ${cont_contigs} -o "${sample_name}_${bins.baseName}_clean.fa"
    """
}