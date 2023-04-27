/*
 * EukCC
*/
process EUKCC {

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    path contigs

    output:
    path "out", emit: contigs

    script:
    """
    bash sed.sh ${contigs}
    """
}

process LINKTABLE {

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    path contigs

    output:
    path "out", emit: contigs

    script:
    """
    bash sed.sh ${contigs}
    """
}