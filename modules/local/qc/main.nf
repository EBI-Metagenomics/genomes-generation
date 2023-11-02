process BUSCO_EUKCC_QC {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    path(eukcc_quality_combined)
    path(busco_list)
    path(genomes_list)
    
    output:
    path("busco_eukcc.qc.csv")
    path("eukcc_final_qc.csv")
    path("busco_final_qc.csv")

    script:
    """
    create_qc_table.py --eukcc_concat ${eukcc_quality_combined} --busco_files ${busco_list} \
    --output "busco_eukcc.qc.csv" --output_eukcc "eukcc_final_qc.csv" --output_busco "busco_final_qc.csv" \
    --genomes_list ${genomes_list}
    """
}
