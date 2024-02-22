
process CONCOCT_CONCOCT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::concoct=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/concoct:1.1.0--py38h7be5676_2':
        'biocontainers/concoct:1.1.0--py38h7be5676_2' }"

    input:
    tuple val(meta), path(coverage_file), path(fasta)

    output:
    tuple val(meta), path("*_args.txt")                     , optional: true,     emit: args_txt
    tuple val(meta), path("*_clustering_gt*.csv"),                                emit: clustering_csv
    tuple val(meta), path("*_log.txt")                      , optional: true,     emit: log_txt
    tuple val(meta), path("*_original_data_gt*.csv")        , optional: true,     emit: original_data_csv
    tuple val(meta), path("*_PCA_components_data_gt*.csv")  , optional: true,     emit: pca_components_csv
    tuple val(meta), path("*_PCA_transformed_data_gt*.csv") , optional: true,     emit: pca_transformed_csv
    path "versions.yml",                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concoct: \$(echo \$(concoct --version 2>&1) | sed 's/concoct //g' )
    END_VERSIONS

    touch ${prefix}_clustering_gt.csv

    set +e
    concoct \\
        $args \\
        --threads ${task.cpus} \\
        --coverage_file ${coverage_file} \\
        --composition_file ${fasta} \\
        -b ${prefix}
    CONCOCT_EXITCODE="\$?"
    set -e

    echo "Exit code \$CONCOCT_EXITCODE"
    if [ "\$CONCOCT_EXITCODE" == "255" ]; then
        echo "CONCOCT exit code \$CONCOCT_EXITCODE"
        if [ ! -e "${prefix}_log.txt" ]; then
            echo "concoct logfile does not exist. Exit"
            echo "Error" >&2
            exit \$CONCOCT_EXITCODE
        else
            echo "concoct logfile exists -> checking for unbinning issue"
            if grep -q "Not enough contigs pass the threshold filter." "${prefix}_log.txt"; then
                echo "Not enough contigs pass the threshold filter"
                exit 0
            else
                echo "Unknown problem. Exit"
                exit \$CONCOCT_EXITCODE
            fi
        fi
    fi
    """
}
