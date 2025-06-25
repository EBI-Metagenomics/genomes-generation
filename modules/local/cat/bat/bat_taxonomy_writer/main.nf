process BAT_TAXONOMY_WRITER {

    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    path bat_names

    output:
    path "human_readable.taxonomy.csv", emit: all_bin2classification_human_readable
    path "all_bin2classification.txt" , emit: all_bin2classification
    path "versions.yml"               , emit: versions

    script:
    """
    bat_taxo_process.py --bat_names ${bat_names} --output "human_readable.taxonomy.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}
