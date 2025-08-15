process GENERATE_INPUT_SAMPLESHEET {

    label 'process_low'

    container "quay.io/microbiome-informatics/ggp_samplesheet_generation:v1"

    secret 'ENA_API_USER'
    secret 'ENA_API_PASSWORD'

    input:
    val(ena_assembly_study_accession)
    val(ena_raw_reads_study_accession)

    output:
    path "*/samplesheet.csv"  , emit: samplesheet
    path "versions.yml"       , emit: versions

    script:
    def scientific_name = params.filter_samplesheet_by_scientific_name ? "--scientific-name ${params.filter_samplesheet_by_scientific_name}" : ""
    def environment_biome = params.filter_samplesheet_by_environment_biome ? "--environment-biome ${params.filter_samplesheet_by_environment_biome}" : ""
    def keep_metat = params.keep_metat_runs_in_samplesheet ? "--keep-metat" : ""
    def args     = task.ext.args ?: ''

    """
    generate_inputs.py \
      --assembly-study ${ena_assembly_study_accession} \
      --raw-reads-study ${ena_raw_reads_study_accession} \
      ${scientific_name} \
      ${environment_biome} \
      ${keep_metat} \
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
