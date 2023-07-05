process QC_BUSCO_EUKCC {

    tag "${name} ${}"

    publishDir(
        path: "${params.outdir}/final_qc_euk",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/biocontainers/biopython:1.75'

    input:
        path eukcc_concoct
	path eukcc_metabat
        path busco_list

    output:
        path "qc.csv"

    script:
    """
    create_qc_table.py --eukcc_c ${eukcc_concoct} --eukcc_m ${eukcc_metabat} --busco_files ${busco_list}
    """
}
