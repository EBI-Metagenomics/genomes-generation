process METAWRAP_BINNING {

    tag "${name}"
    container 'quay.io/microbiome-informatics/metawrap:latest'

    publishDir(
        path: "${params.outdir}/${name}/binning/",
        mode: 'copy',
        failOnError: true
    )

    input:
    tuple val(name), path(contigs)
    path input_reads
    val binner

    output:
    tuple val(name), path("output_${binner}"), emit: binning
    tuple val(name), path("binning/work_files/metabat_depth.txt"), emit: metabat_depth_for_coverage, optional: true

    script:
    reads = input_reads.collect()
    def args = " --${binner} ";
    if ( input_reads.size() == 1 ) {
        args += "--single-end ${input_reads}"
    }
    if ( input_reads.size() == 2 ) {
        args += "${input_reads[0]} ${input_reads[1]}"
    }
    """
    echo "Running binning"
    metawrap binning -t ${task.cpus} -m 30 -l 2500 -a ${contigs} -o binning ${args}

    mkdir -p output_${binner}
    cd binning/*_bins
    for f in \$(ls . | grep ".fa") ;
    do
        mv "\$f" "../../output_${binner}/${name}_${binner}_\$f" ;
    done
    """

    stub:
    """
    mkdir -p output_${binner} binning/work_files
    touch binning/work_files/metabat_depth.txt
    """
}

process BIN_REFINEMENT {

    //container 'quay.io/microbiome-informatics/metawrap:latest'
    tag "${name}"

    publishDir(
        path: "${params.outdir}/binning/",
        mode: 'copy',
        failOnError: true
    )

    input:
    tuple val(name), path(binning_metabat2), path(binning_concoct), path(binning_maxbin2)

    output:
    tuple val(name), path("bin_refinement/metawrap_*bins/*"), emit: bin_ref_bins

    script:
    """
    echo "Running bin_refinement"

    metawrap bin_refinement -t ${task.cpus} -o bin_refinement \
    -A ${binning_metabat2} -B ${binning_concoct} -C ${binning_maxbin2} \
    -c 50 -x 5 -m ${task.memory}
    """

    stub:
    """
    mkdir -p bin_refinement/metawrap_50_5_bins
    touch bin_refinement/metawrap_50_5_bins/${name}.bin.0.fa bin_refinement/metawrap_50_5_bins/${name}.bin.1.fa
    """
}

