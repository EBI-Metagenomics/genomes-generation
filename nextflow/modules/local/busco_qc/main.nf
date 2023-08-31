process BUSCO_QC {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
        path(eukcc_quality_combined)
        path(busco_list)

    output:
        path("qc.csv")

    script:
    """
    create_qc_table.py --eukcc_concat ${eukcc_quality_combined} --busco_files ${busco_list} --output "qc.csv"
    """
}
