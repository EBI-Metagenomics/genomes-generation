/*
 * EukCC
*/
process EUKCC {

    publishDir(
        path: "${params.outdir}/eukcc/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    val binner
    tuple val(name), path(links)
    path eukcc_db
    tuple val(name), path(bindir)

    output:
    tuple val(name), path("*_merged_bins"), emit: eukcc_results
    tuple val(name), path("eukcc.csv"), emit: eukcc_csv

    script:
    """
    eukcc folder \
        --improve_percent 10 \
        --n_combine 1 \
        --threads ${task.cpus} \
        --improve_ratio  5 \
        --links ${links} \
        --min_links 100 \
        --suffix .fa \
        --db ${eukcc_db} \
        --out ${binner}_${name}_merged_bins \
        --prefix "${binner}_${name}_merged." \
        ${bindir}
    """
}

process LINKTABLE {

    publishDir(
        path: "${params.outdir}/eukcc/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/microbiome-informatics/eukrecover.python3base:v1'

    input:
    tuple val(name), path(bindir)
    tuple val(name), path(bam)

    output:
    tuple val(name), path("*.links.csv"), emit: links_table

    script:
    """
    binlinks.py  --ANI 99 --within 1500 --out ${name}.links.csv ${bindir} ${bam}
    """
}