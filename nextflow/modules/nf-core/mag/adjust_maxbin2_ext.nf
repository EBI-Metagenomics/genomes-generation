process ADJUST_MAXBIN2_EXT {
    tag "${meta.id}"
    label 'process_low'

    // Using container from multiqc since it'll be included anyway
    //conda "bioconda::multiqc=1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("${meta.id}_maxbin_bins"), emit: renamed_bins

    script:
    """
    mkdir -p "${meta.id}_maxbin_bins"

    if [ -n "${bins}" ]
    then
        for file in ${bins}; do
            [[ \${file} =~ (.*).fasta ]];
            bin="\${BASH_REMATCH[1]}"
            mv \${file} \${bin}.fa
        done
    fi
    mv *.fa "${meta.id}_maxbin_bins"
    """
}
