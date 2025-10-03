process CONCOCT_EXTRACTFASTABINS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::concoct=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/concoct:1.1.0--py38h7be5676_2':
        'quay.io/biocontainers/concoct:1.1.0--py38h7be5676_2' }"

    input:
    tuple val(meta), path(original_fasta), path(csv)

    output:
    tuple val(meta), path("${meta.id}_concoct_bins")     , emit: fastas_dir, optional: true
    tuple val(meta), path("${meta.id}_concoct_bins/*.fa"), emit: fastas, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    extract_fasta_bins.py \\
        $args \\
        $original_fasta \\
        $csv \\
        --output_path ${prefix}

    export BINS=\$(ls ${prefix} | wc -l)
    if [ \$BINS -eq 0 ]; then
        echo "Result folder is empty"
    else
        ## Add prefix to each file to disambiguate one sample's 1.fa, 2.fa from sample2
        ## renames 1.fa to accession_1.fa
        for i in ${prefix}/*.fa; do
            mv \${i} \${i/\\///${prefix}_}
        done

        ## Create final output folder
        ## This helps to avoid output folder being empty
        mkdir -p ${meta.id}_concoct_bins

        ## Move files to the output folder
        for i in ${prefix}/*.fa; do
            mv \${i} ${meta.id}_concoct_bins/
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concoct: \$(echo \$(concoct --version 2>&1) | sed 's/concoct //g' )
    END_VERSIONS
    """
}