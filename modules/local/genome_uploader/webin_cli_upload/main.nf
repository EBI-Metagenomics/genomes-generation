process WEBIN_CLI_UPLOAD {

    label 'process_low'
    tag "${id}"
    stageInMode 'copy'
    container "quay.io/biocontainers/ena-webin-cli:8.2.0--hdfd78af_0"

    input:
    tuple val(id), path(mag), path(manifest)

    output:
    tuple val(id), path("*manifest.report") , emit: webin_report
    path "versions.yml"                     , emit: versions

    script:
    """
    # change FASTA path in manifest to current workdir
    export MAG_FULL_PATH=\$(readlink -f ${mag})
    sed 's|^FASTA\t.*|FASTA\t'"\${MAG_FULL_PATH}"'|g' ${manifest} > ${id}_updated_manifest.manifest

    ena-webin-cli \
      -context=genome \
      -manifest=${id}_updated_manifest.manifest \
      -userName='$params.webin_account' \
      -password='$params.webin_password' \
      -submit

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ena-webin-cli: \$(ena-webin-cli -version 2>&1 )
    END_VERSIONS
    """
}
