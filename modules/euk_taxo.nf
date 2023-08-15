process EUK_TAXONOMY_WRITER {

    tag "${bat_names}"

    publishDir(
        path: "${params.outdir}/euk_taxonomy",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/biocontainers/biopython:1.75'

    input:
        path bat_names

    output:
        path "euk_taxonomy.csv"
	path "all_bin2classification.txt"

    script:
    """
    bat_taxo_process.py --bat_names ${bat_names}
    """
}
