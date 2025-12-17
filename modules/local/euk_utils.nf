process CONCATENATE_QUALITY_FILES {
    tag "${meta.id}"

    input:
    tuple val(meta), path(input_files)
    val output_name

    output:
    tuple val(meta), path("${meta.id}.${output_name}"), emit: concatenated_result

    script:
    """
    echo "bin\tcompleteness\tcontamination\tncbi_lng" > "${meta.id}.${output_name}"
    for i in ${input_files}; do
        tail -n +2 \${i} > help_file
        if [ -s help_file ]; then
            grep -v "completeness" \${i} >> "${meta.id}.${output_name}"
        fi
    done
    rm help_file
    """
}

process MODIFY_QUALITY_FILE {
    tag "${output_name}"

    input:
    path(quality_table_csv)
    val output_name

    output:
    path("${output_name}"), emit: modified_result

    script:
    """
    echo "genome,completeness,contamination" > ${output_name}
    grep -v "completeness" ${quality_table_csv} | cut -f1-3 | tr '\t' ',' >> ${output_name} || true
    """
}