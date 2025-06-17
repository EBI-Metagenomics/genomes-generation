/*
 * change . to _ from contigs
*/
process MODIFY_CONTIGS {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip:416534c89faa8870':
        'community.wave.seqera.io/library/pip_base:fa76cdb8dd2b80d6' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_underscore.fasta.gz"), emit: underscore_contigs

    script:
    """
    zcat ${contigs} | sed 's/\\./\\_/' | gzip > ${meta.id}_underscore.fasta.gz
    """
}