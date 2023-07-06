process QC_BUSCO_EUKCC {

    publishDir(
        path: "${params.outdir}/final_qc_euk",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/biocontainers/biopython:1.75'

    input:
        path eukcc_combined
        path busco_list

    output:
        path "qc.csv"

    script:
    """
    create_qc_table.py --eukcc_concat ${eukcc_combined} --busco_files ${busco_list}
    """
}
