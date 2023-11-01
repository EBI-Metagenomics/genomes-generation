process RENAME_MAXBIN {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::maxbin2=2.2.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/maxbin2:2.2.7--he1b5a44_2' :
        'quay.io/biocontainers/maxbin2:2.2.7--he1b5a44_2' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("${meta.id}_maxbin_bins"), emit: renamed_bins

    errorStrategy {
        task.exitStatus {
            exitVal ->
                // Retry on non-zero exit codes
                return exitVal != 0 ? ErrorAction.RETRY : ErrorAction.FINISH
        }
        maxRetries 3  // Set the maximum number of retries
        sleep 10       // Set the delay between retries in seconds
    }

    script:
    """
    version=\$( run_MaxBin.pl -v | head -n 1 | sed 's/MaxBin //' )

    rename_maxbin.py -o ${meta.id}_maxbin_bins -v \${version} -a ${meta.id} --bins ${bins}
    """
}