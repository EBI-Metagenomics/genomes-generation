process EUK_TAXONOMY {
    tag "${bin}"

    container 'quay.io/microbiome-informatics/cat:5.2.3'

    publishDir(
        path: "${params.outdir}/bat_output",
        mode: 'copy',
        failOnError: true
    )

    input:
    path bin
    path cat_db
    path taxonomy_db

    output:
    path '*.BAT_run.bin2classification.names.txt', emit: bat_names

    script:
    """
    echo "[MAG euk taxonomy] Analysing bins"
    CAT bin -b ${bin} -d ${cat_db} -t ${taxonomy_db} -o ${bin.baseName}.BAT_run --force

    echo "[MAG euk taxonomy] Adding taxonomy names"
    CAT add_names -i ${bin.baseName}.BAT_run.bin2classification.txt -o ${bin.baseName}.BAT_run.bin2classification.names.txt -t ${taxonomy_db} --only_official

    """
}
