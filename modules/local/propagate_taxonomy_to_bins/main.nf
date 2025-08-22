process PROPAGATE_TAXONOMY_TO_BINS {
    tag "${meta.id}"
    label 'process_single'

    container ''

    input:
    tuple val(meta), path(drep_clustering_csv)
    path(taxonomy_tsv)

    output:
    tuple val(meta), path("ncbi_taxonomy.txt"), emit: ncbi_taxonomy
    path "versions.yml"                       , emit: versions


    script:
    """
    propagate_taxonomy_to_bins.py \
        --cdb ${drep_clustering_csv} \
        --taxonomy ${taxonomy_tsv} \
        --output ncbi_taxonomy_propagated.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """

    stub:
    """
    touch ncbi_taxonomy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}