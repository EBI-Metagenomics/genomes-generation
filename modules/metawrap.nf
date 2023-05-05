process METAWRAP_BINNING {

    container 'quay.io/microbiome-informatics/metawrap:latest'

    publishDir(
        path: "${params.outdir}/binning/",
        mode: 'copy',
        failOnError: true
    )

    input:
    val mode
    val name
    path contigs
    path input_reads

    output:
    path "binning/metabat2_bins", emit: binning_metabat2
    path "binning/concoct_bins", emit: binning_concoct
    path "binning/maxbin2_bins", emit: binning_maxbin2

    script:
    def args = "";
    if ( mode == "single" ) {
        args = "--single-end ${input_reads}"
    }
    if ( mode == "paired" ) {
        args = "${input_reads[0]} ${input_reads[1]}"
    }
    """
    echo "Running binning"
    metawrap binning -t 8 -m 10 -l 2500 --metabat2 --concoct --maxbin2 -a ${contigs} -o binning ${args}
    """
}

process BIN_REFINEMENT {

    container 'quay.io/microbiome-informatics/metawrap:latest'

    publishDir(
        path: "${params.outdir}/binning/",
        mode: 'copy',
        failOnError: true
    )

    input:
    path binning_metabat2
    path binning_concoct
    path binning_maxbin2

    output:
    path "bin_refinement", emit: bin_ref

    script:
    """
    echo "Running bin_refinement"

    metawrap bin_refinement -t ${task.cpus} -o bin_refinement \
    -A ${binning_metabat2} -B ${binning_concoct} -C ${binning_maxbin2} \
    -c 50 -x 5 -m ${task.memory}
    """
}

