/*
 * EukCC
*/
process EUKCC {

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    val name
    path links
    path eukcc_db
    path bindir

    output:
    path "output", emit: eukcc_results

    script:
    """
    eukcc \
        --improve_percent 10 \
        --n_combine 1 \
        --threads ${task.cpus} \
        --improve_ratio  5 \
        --links ${links} \
        --min_links 100 \
        --suffix .fa \
        --db ${eukcc_db} \
        --out output \
        --prefix "${name}_merged." \
        ${bindir}
    """
}

process LINKTABLE {

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    val name
    path bindir
    path bam

    output:
    path "out", emit: contigs

    script:
    """
    python3 binlinks.py  --ANI 99 --within 1500 --out ${name}.links.csv ${bindir} ${bam}
    """
}