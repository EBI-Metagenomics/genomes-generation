process BUSCO_EUKCC_QC {

    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    path eukcc_quality_combined
    path busco_list
    path genomes_list
    
    output:
    path "combined_busco_eukcc.qc.csv", emit: busco_eukcc
    path "eukcc_final_qc.csv", emit: eukcc_final_qc
    path "busco_final_qc.csv", emit: busco_final_qc
    path "versions.yml"      , emit: versions

    script:
    """
    create_qc_table.py \
        --eukcc_concat ${eukcc_quality_combined} \
        --busco_files ${busco_list} \
        --output "combined_busco_eukcc.qc.csv" \
        --output_eukcc "eukcc_final_qc.csv" \
        --output_busco "busco_final_qc.csv" \
        --genomes_list ${genomes_list}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
