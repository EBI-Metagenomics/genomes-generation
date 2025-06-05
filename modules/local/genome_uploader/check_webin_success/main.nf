process CHECK_WEBIN_SUCCESS {
    input:
    tuple val(id), path(report_file)

    output:
    tuple val(id), env(success_status), emit: results

    script:
    """
    if grep -q "submission has been completed successfully" ${report_file}; then
        export success_status="true"
    else
        export success_status="false"
    fi
    """
}