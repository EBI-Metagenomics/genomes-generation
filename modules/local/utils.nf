/*
    ~~~~~~~~~~~~~~~~~~
     uncompress
    ~~~~~~~~~~~~~~~~~~
*/
process GUNZIP {
    input:
    path(compressed_file)

    output:
    path("out/*"), emit: uncompressed

    script:
    """
    mkdir out
    cp ${compressed_file} out
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

    tag "${meta.id}"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${contigs}"), emit: underscore_contigs

    script:
    """
    sed -i 's/\\./\\_/' ${contigs}
    """
}

process ERR_TO_ERZ_SED {

    tag "${meta.id}"

    container 'quay.io/biocontainers/pigz:2.3.4'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*changed*.*.gz"), emit: modified_reads
    path "versions.yml"                       , emit: versions

    script:
    input_ch = reads.collect()
    if (input_ch.size() == 1 ) {
        """
        zcat "${input_ch[0]}" | sed 's/${meta.id}/${meta.erz}/g' | pigz > ${meta.id}_changed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pigz: \$(pigz --version 2>&1 | sed 's/pigz //g')
        END_VERSIONS
        """
    }
    else if (input_ch.size() == 2 ) {
        """
        zcat "${input_ch[0]}" | sed 's/${meta.id}/${meta.erz}/g' | pigz > ${meta.id}_changed_1.fastq.gz
        zcat "${input_ch[1]}" | sed 's/${meta.id}/${meta.erz}/g' | pigz > ${meta.id}_changed_2.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pigz: \$(pigz --version 2>&1 | sed 's/pigz //g')
        END_VERSIONS
        """
    }
}

/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_UNDERSCORE_TO_DOT {

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

