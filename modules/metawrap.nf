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
    def reads = "";
    def file = "";
    if ( mode == "single" ) {
        args = "--single-end"
        reads = "${input_reads}"
        file = "${input_reads}"
    }
    if ( mode == "paired" ) {
        reads = "${input_reads[0]} ${input_reads[1]}"
        file = "${input_reads[0]}"
    }
    """
    if (file ${file} | grep -q compressed ) ; then
        echo "gunzip"
        gunzip ${reads}
    fi
    echo "Running binning"
    metawrap binning ${args} -t 8 -m 10 -l 2500 --metabat2 --concoct --maxbin2 -a ${contigs} -o binning ${name}*.f*q
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
    path "binning", emit: bin_ref

    script:
    def args = "";
    print("${binning_metabat2}")
    if (file("${binning_metabat2}").list() == []) {
        args += " -A ${binning_metabat2}"
    }

    if (file("${binning_concoct}").list() == []) {
        args += " -B ${binning_concoct}"
    }

    if (file("${binning_maxbin2}").list() ) {
        args += " -C ${binning_maxbin2}"
    }
    """
    echo "Running bin_refinement"

    metawrap bin_refinement -t ${task.cpus} -o bin_refinement \
    "${args}" \
    -c 50 -x 5 -m ${task.memory}
    """
}

