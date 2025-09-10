process GTDBTK_TO_NCBI_TAXONOMY {

    label 'process_long'
    container 'quay.io/biocontainers/gtdbtk:2.4.1--pyhdfd78af_1'

    containerOptions "${ workflow.containerEngine == 'singularity' ?
        "--bind ${gtdbtk_refdata}:/opt/gtdbtk_refdata":
        "-v ${gtdbtk_refdata}:/opt/gtdbtk_refdata" }"

    input:
    path(gtdbtk_results)
    path(gtdbtk_refdata)

    output:
    path "ncbi_taxonomy.txt", emit: ncbi_taxonomy
    path "versions.yml"     , emit: versions


    script:
    """
    export GTDBTK_DATA_PATH=/opt/gtdbtk_refdata

    echo "Create NCBI taxonomy"
    gtdb_to_ncbi_majority_vote.py \
        --gtdbtk_output_dir ${gtdbtk_results} \
        --output_file ncbi_taxonomy.txt \
        --ar53_metadata_file ${gtdbtk_refdata}/ar53_metadata_r???.tsv \
        --bac120_metadata_file ${gtdbtk_refdata}/bac120_metadata_r???.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb_to_ncbi_majority_vote.py: \$(echo \$(gtdb_to_ncbi_majority_vote.py -h 2>&1) | head -n 1 | sed "s/gtdb_to_ncbi_majority_vote.py v//; s/:.*//")
    END_VERSIONS
    """

    stub:
    """
    export GTDBTK_DATA_PATH=/opt/gtdbtk_refdata

    touch ncbi_taxonomy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdb_to_ncbi_majority_vote.py: \$(echo \$(gtdb_to_ncbi_majority_vote.py -h 2>&1) | head -n 1 | sed "s/gtdb_to_ncbi_majority_vote.py v//; s/:.*//")
    END_VERSIONS
    """
}