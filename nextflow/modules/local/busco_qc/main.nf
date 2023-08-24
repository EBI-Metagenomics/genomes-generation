process BUSCO_QC {

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
        tuple val(meta), path(eukcc_quality_combined)
        path(busco_list)

    output:
        tuple val(meta), path("*qc.csv")

    script:
    """
    create_qc_table.py --eukcc_concat ${eukcc_quality_combined} --busco_files ${busco_list} --output "${meta.id}.qc.csv"
    """
}
