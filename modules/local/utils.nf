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

    //publishDir(
    //    path: "${params.outdir}/${output_folder}",
    //    mode: 'copy',
    //    failOnError: true
    //)

    input:
    tuple val(meta), file(file_to_compress)
    //val output_folder

    output:
    path("*.gz"), emit: compressed

    script:
    """
    pigz ${file_to_compress}
    """
}

/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_DOT_TO_UNDERSCORE_CONTIGS {

    tag "${meta.id}"
    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${contigs}"), emit: underscore_contigs

    script:
    """
    sed -i 's/\\./\\_/' ${contigs}
    """
}

process ERZ_TO_ERR {

    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(reads)
    path rename_file

    output:
    tuple val(meta), path("changed*.fastq.gz"), emit: modified_reads

    script:
    """
    echo ${meta.id}
    grep "${meta.id}" ${rename_file} > help_file
    export from_accession=\$(cat help_file | cut -f1)
    export to_accession=\$(cat help_file | cut -f2)
    echo "\${from_accession} --> \${to_accession}"

    gunzip ${reads}

    change_reads.py --reads *.fastq -f \${from_accession} -t \${to_accession} --change_dots_to_underscores

    gzip changed*.fastq
    """
}

/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_UNDERSCORE_TO_DOT {

    tag "${name} ${contigs}"

    input:
    path(contigs)

    output:
    path("${contigs}"), emit: return_files

    script:
    """
    sed -i 's/\\_/\\./' ${contigs}
    """
}

