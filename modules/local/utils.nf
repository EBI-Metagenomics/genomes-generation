/*
    ~~~~~~~~~~~~~~~~~~
     uncompress
    ~~~~~~~~~~~~~~~~~~
*/
process GUNZIP {

    label 'process_low'
    tag "$meta.id"

    input:
    tuple val(meta), path(compressed_file)

    output:
    tuple val(meta), path("out/*"), emit: uncompressed

    script:
    """
    mkdir out
    cp ${compressed_file} "out/${meta.id}.fasta.gz"
    cd out
    gunzip *
    """
}

/*
    ~~~~~~~~~~~~~~~~~~
     compress
    ~~~~~~~~~~~~~~~~~~
*/
process GZIP {

    label 'process_low'
    stageInMode 'copy'
    tag "${file_to_compress}"

    input:
    file(file_to_compress)

    output:
    path("*.gz"), emit: compressed
    path "versions.yml"                       , emit: versions

    script:
    """
    pigz ${file_to_compress}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """
}


/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_UNDERSCORE_TO_DOT {

    label 'process_low'
    tag "${contigs}"

    input:
    path(contigs)

    output:
    path("${contigs}"), emit: return_files

    script:
    """
    sed -i 's/\\_/\\./' ${contigs}
    """
}


process FINALIZE_LOGGING {

    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:0.24.1':
        'quay.io/biocontainers/pandas:0.24.1' }"

    input:
    path(logging_file)
    val(output)

    output:
    path("${output}"), emit: structured_logging

    script:
    """
    logging_stats.py -i ${logging_file} -o ${output}
    """
}
