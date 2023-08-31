process BAT_TAXONOMY_WRITER {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
        path(bat_names)

    output:
        path("*euk_taxonomy.csv")
        path("all_bin2classification.txt")

    script:
    """
    bat_taxo_process.py --bat_names ${bat_names} --output "${bat_names.baseName}.euk_taxonomy.csv"
    """
}
