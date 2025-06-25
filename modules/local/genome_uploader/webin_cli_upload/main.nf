process WEBIN_CLI_UPLOAD {

    label 'process_low'
    tag "${id}"
    stageInMode 'copy'
    container "quay.io/biocontainers/ena-webin-cli:8.2.0--hdfd78af_0"

    input:
    secret 'ENA_API_USER'
    secret 'ENA_API_PASSWORD'
    tuple val(id), path(mag), path(manifest)

    output:
    tuple val(id), path("webin-cli.report") , emit: webin_report
    tuple val(id), env('SUCCESS_STATUS')    , emit: upload_status
    path "versions.yml"                     , emit: versions

    script:
    """
    # change FASTA path in manifest to current workdir
    export MAG_FULL_PATH=\$(readlink -f ${mag})
    sed 's|^FASTA\t.*|FASTA\t'"\${MAG_FULL_PATH}"'|g' ${manifest} > ${id}_updated_manifest.manifest

    ena-webin-cli \
      -context=genome \
      -manifest=${id}_updated_manifest.manifest \
      -userName='\$ENA_API_USER' \
      -password='\$ENA_API_PASSWORD' \
      -submit

    if grep -q "submission has been completed successfully" webin-cli.report; then
        export SUCCESS_STATUS="true"
    else
        export SUCCESS_STATUS="false"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ena-webin-cli: \$(ena-webin-cli -version 2>&1 )
    END_VERSIONS
    """
}
