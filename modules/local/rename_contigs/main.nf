/*
 * change headers to short names for contigs fasta
*/
process RENAME_CONTIGS {

    label 'process_low'
    tag "${meta.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_renamed.fasta.gz"), emit: contigs_renamed
    path "versions.yml",                                  emit: versions

    script:
    """
    rename_contigs \\
       --input ${contigs} \\
       --output ${meta.id}_renamed.fasta \\
       --prefix ${meta.assembly_accession}
    gzip ${meta.id}_renamed.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}