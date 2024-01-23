process GTDBTK_TO_NCBI_TAXONOMY {

    label 'process_long'
    container 'quay.io/microbiome-informatics/gtdb-tk:2.3.0.1'

    containerOptions "${ workflow.containerEngine == 'singularity' ?
        "--bind ${gtdbtk_refdata}:/opt/gtdbtk_refdata":
        "-v ${gtdbtk_refdata}:/opt/gtdbtk_refdata" }"

    input:
    path(gtdbtk_results)
    path(gtdbtk_refdata)

    output:
    path "ncbi_taxonomy.txt", emit: ncbi_taxonomy
    path "versions.yml"         , emit: versions


    script:
    """
    export GTDBTK_DATA_PATH=/opt/gtdbtk_refdata

    echo "Create NCBI taxonomy"
    gtdb_to_ncbi_majority_vote.py \
        --gtdbtk_output_dir ${gtdbtk_results} \
        --output_file ncbi_taxonomy.txt \
        --ar53_metadata_file ${gtdbtk_refdata}/ar53_metadata_r214.tsv \
        --bac120_metadata_file ${gtdbtk_refdata}/bac120_metadata_r214.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb_to_ncbi_majority_vote.py: \$(echo \$(gtdb_to_ncbi_majority_vote.py --version 2>&1) | grep -e "INFO: GTDB to NCBI majority vote v" | sed "s/INFO: GTDB to NCBI majority vote v//")
    END_VERSIONS
    """

    stub:
    """
    export GTDBTK_DATA_PATH=/opt/gtdbtk_refdata

    touch ncbi_taxonomy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb_to_ncbi_majority_vote.py: \$(echo \$(gtdb_to_ncbi_majority_vote.py --version 2>&1) | grep -e "INFO: GTDB to NCBI majority vote v" | sed "s/INFO: GTDB to NCBI majority vote v//")
    END_VERSIONS
    """
}
