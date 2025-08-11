process CREATE_MANIFESTS_FOR_UPLOAD {

    label 'process_low'

    container "community.wave.seqera.io/library/genome-uploader:2.4.0--c32667b153eea137"

    input:
    secret 'WEBIN_ACCOUNT'
    secret 'WEBIN_PASSWORD'
    path(table_for_upload)
    path(mags)

    output:
    path "results/MAG_upload/manifests*/*.manifest", emit: manifests
    path "results/MAG_upload/ENA_backup.json"      , emit: ena_upload_backup_json
    path "results/MAG_upload/genome_samples.xml"   , emit: upload_genome_samples
    path "results/MAG_upload/registered_MAGs*.tsv" , emit: upload_registered_mags
    path "results/MAG_upload/submission.xml"       , emit: upload_submission_xml
    //path "versions.yml"       , emit: versions

    script:
    def mags_arg = params.upload_mags ? "--mags" : ""
    def bins_arg = params.upload_bins ? "--bins" : ""
    def tpa      = params.upload_tpa  ? "--tpa"  : ""
    def force    = params.upload_force  ? "--force"  : ""
    def mode     = (!params.test_upload) ? "--live" : ""
    def args     = task.ext.args ?: ''

    """
    export ENA_WEBIN=\$WEBIN_ACCOUNT
    export ENA_WEBIN_PASSWORD=\$WEBIN_PASSWORD

    genome_upload \
      -u $params.ena_assembly_study_accession \
      --genome_info ${table_for_upload} \
      --centre_name $params.centre_name \
      ${mags_arg} \
      ${bins_arg} \
      ${tpa} \
      ${force} \
      ${mode} \
      --out results \
      ${args}
    """
}
