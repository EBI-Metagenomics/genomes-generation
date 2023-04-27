process BINNING {

    container 'quay.io/microbiome-informatics/metawrap:latest'

    publishDir(
        path: "${params.outdir}/",
        mode: 'copy',
        failOnError: true
    )

    input:
    val mode
    path contigs
    val fastq_single
    path fastq_pf
    path fastq_pr

    output:
    path "binning/metabat2_bins", emit: binning_metabat2
    path "binning/concoct_bins", emit: binning_concoct
    path "binning/maxbin2_bins", emit: binning_maxbin2

    script:
    def args = "";
    def reads = "";
    if ( mode == "single" ) {
        args = "--single-end"
        reads = "${fastq_single}"
    }
    if ( mode == "paired" ) {
        reads = "${fastq_pf} ${fastq_pr}"
    }
    """
    echo "Running binning"
    metawrap binning ${args} -t 8 -m 10 -l 2500 --metabat2 --concoct --maxbin2 -a ${contigs} -o binning ${reads}
    """
}


process BIN_REFINEMENT {

    container 'quay.io/microbiome-informatics/metawrap:latest'

    publishDir(
        path: "${params.outdir}/",
        mode: 'copy',
        failOnError: true
    )

    input:
    val mode_analysis
    path binning_metabat2
    path binning_concoct
    path binning_maxbin2

    output:
    path "binning", emit: bin_ref

    script:
    def args = '-A ${binning_metabat2} -B ${binning_concoct}';
    if (mode_analysis == 'prok') {
        args += '-C ${binning_maxbin2}'
    }
    """
    echo "Running bin_refinement"
    metawrap bin_refinement -t 4 -o bin_refinement ${args} -c 50 -x 5
    """
}