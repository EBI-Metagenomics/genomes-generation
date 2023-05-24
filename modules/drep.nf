/*
 * dRep, this workflow dereplicates a set of genomes.
*/
process DREP {

    publishDir(
        path: "${params.outdir}/drep",
        mode: 'copy',
    )

    container 'quay.io/biocontainers/drep:3.2.2--pyhdfd78af_0'

    input:
    tuple val(name), path(genomes_directory)
    tuple val(name), path(quality_csv)
    val nc
    val pa
    val sa
    val comp
    val con

    output:
    path "drep_output/dereplicated_genomes/*", emit: dereplicated_genomes

    script:
    """
    dRep dereplicate -g ${genomes_directory}/*.fa \
    -p ${task.cpus} \
    -pa ${pa} \
    -sa ${sa} \
    -nc ${nc} \
    -cm larger \
    -comp ${comp} \
    -con ${con} \
    --genomeInfo ${quality_csv} \
    drep_output
    """
}
