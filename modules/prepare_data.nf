/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_DOT_TO_UNDERSCORE {

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    publishDir(
        path: "${params.outdir}/prepare_data",
        mode: 'copy',
        failOnError: true
    )

    input:
    tuple val(accession), path(contigs)

    output:
    tuple val(accession), path("${contigs}"), emit: return_contigs

    script:
    """
    sed -i 's/\\./\\_/' ${contigs}
    """
}

/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_UNDERSCORE_TO_DOT {

    input:
    path contigs

    output:
    path "${contigs.baseName}", emit: return_contigs

    script:
    """
    sed -i 's/\\_/\\./' ${contigs}
    """
}