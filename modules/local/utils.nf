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
process CHANGE_DOT_TO_UNDERSCORE_CONTIGS {

    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_underscore.fasta.gz"), emit: underscore_contigs

    script:
    """
    zcat ${contigs} | sed 's/\\./\\_/' | gzip > ${meta.id}_underscore.fasta.gz
    """
}

process ERR_TO_ERZ {

    label 'process_medium'
    tag "${meta.id}"

    container 'quay.io/biocontainers/seqkit:2.7.0--h9ee0642_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*sanitized*.*.gz"), emit: modified_reads
    path "versions.yml"                       , emit: versions

    script:
    input_ch = reads.collect()
    if (meta.single_end ) {
        """
        seqkit sana ${input_ch[0]} | seqkit replace -p ${meta.id} -r ${meta.erz} -o ${meta.id}_sanitized.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version 2>&1 | sed 's/seqkit //g')
        END_VERSIONS
        """
    }
    else {
        """
        echo "read1"
        seqkit sana ${input_ch[0]} | seqkit replace -p ${meta.id} -r ${meta.erz} -o ${meta.id}_sanitized_1.fastq.gz

        echo "read2"
        seqkit sana ${input_ch[1]} | seqkit replace -p ${meta.id} -r ${meta.erz} -o ${meta.id}_sanitized_2.fastq.gz

        echo "sync read1 and read2"
        seqkit pair -1 ${meta.id}_sanitized_1.fastq.gz -2 ${meta.id}_sanitized_2.fastq.gz -O result

        mv result/${meta.id}_sanitized_1.fastq.gz ${meta.id}_sanitized_1.fastq.gz
        mv result/${meta.id}_sanitized_2.fastq.gz ${meta.id}_sanitized_2.fastq.gz
        rm -rf result

        echo "Done"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version 2>&1 | sed 's/seqkit //g')
        END_VERSIONS
        """
    }
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
