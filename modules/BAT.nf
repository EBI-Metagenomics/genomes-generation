process EUK_TAXONOMY {
    tag "${bin}"

    container 'quay.io/microbiome-informatics/cat:5.2.3'

    publishDir(
        path: "${params.outdir}/euk_taxonomy",
        mode: 'copy',
        failOnError: true
    )

    input:
    path bin
    path cat_db
    path taxonomy_db

    output:
    path '*.summary.txt'
    path '*.BAT_run.ORF2LCA.names.txt'

    script:
    """
    echo "[MAG euk taxonomy] Analysing bins"
    CAT bin -b ${bin} -d ${cat_db} -t ${taxonomy_db} -o ${bin.baseName}.BAT_run

    echo "[MAG euk taxonomy] Adding taxonomy names"
    CAT add_names -i ${bin.baseName}.BAT_run.ORF2LCA.txt -o ${bin.baseName}.BAT_run.ORF2LCA.names.txt -t ${taxonomy_db} --only_official

    echo "[MAG euk taxonomy] Summarizing output"
    CAT summarise -i ${bin.baseName}.BAT_run.bin2classification.txt -o ${bin.baseName}.summary.txt

    """
}
