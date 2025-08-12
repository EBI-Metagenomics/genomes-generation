process WEBIN_CLI_UPLOAD {

    label 'process_low'
    tag "${id}"
    stageInMode 'copy'
    container "quay.io/biocontainers/ena-webin-cli:8.2.0--hdfd78af_0"
    errorStrategy = { task.exitStatus in ((130..145) + 104 + 1) ? 'retry' : 'finish' }

    input:
    secret 'WEBIN_ACCOUNT'
    secret 'WEBIN_PASSWORD'
    tuple val(id), path(mag), path(manifest)

    output:
    tuple val(id), path("*webin-cli.report") , emit: webin_report
    tuple val(id), env('STATUS')             , emit: upload_status
    path "versions.yml"                      , emit: versions

    script:

    def mode     = params.test_upload ? "-test" : ""

    """
    # change FASTA path in manifest to current workdir
    export MAG_FULL_PATH=\$(readlink -f ${mag})
    sed 's|^FASTA\t.*|FASTA\t'"\${MAG_FULL_PATH}"'|g' ${manifest} > ${id}_updated_manifest.manifest

    ena-webin-cli \
      -context=genome \
      -manifest=${id}_updated_manifest.manifest \
      -userName='\$WEBIN_ACCOUNT' \
      -password='\$WEBIN_PASSWORD' \
      -submit \
      ${mode}

    mv webin-cli.report "${id}_webin-cli.report"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ena-webin-cli: \$(ena-webin-cli -version 2>&1 )
    END_VERSIONS

    # status check
    if grep -q "submission has been completed successfully" "${id}_webin-cli.report"; then
        # first time submission completed successfully
        export STATUS="success"
        true
    elif grep -q "object being added already exists in the submission account with accession" "${id}_webin-cli.report"; then
        # there was attempt to re-submit already submitted genome
        export STATUS="success"
        true
    else
        export STATUS="failed"
        false
    fi
    """
}
