process BAT {
    tag "${bin}"

    container 'quay.io/biocontainers/cat:5.2.3--hdfd78af_1'

    input:
    path bin
    path cat_db_folder
    path cat_taxonomy_db

    output:
    path '*.BAT_run.bin2classification.names.txt', emit: bat_names
    path "versions.yml"                          , emit: versions

    script:
    """
    echo "[MAG euk taxonomy] Analysing bins"
    CAT bin -b ${bin} \
      -d ${cat_db_folder} \
      -t ${cat_taxonomy_db} \
      -o ${bin.baseName}.BAT_run \
      --force

    echo "[MAG euk taxonomy] Adding taxonomy names"
    CAT add_names -i ${bin.baseName}.BAT_run.bin2classification.txt \
      -o ${bin.baseName}.BAT_run.bin2classification.names.txt \
      -t ${cat_taxonomy_db} \
      --only_official

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CAT: \$(CAT --version | sed "s/CAT v//; s/(.*//")
    END_VERSIONS
    """
}
