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
    path genomes_directory
    path checkm_csv
    val nc

    output:
    path "drep_output/dereplicated_genomes", emit: dereplicated_genomes

    script:
    """
    dRep dereplicate -g ${genomes_directory}/*.fa \
    -p ${task.cpus} \
    -pa 0.9 \
    -sa 0.95 \
    -nc ${nc} \
    -cm larger \
    -comp 50 \
    -con 5 \
    --genomeInfo ${checkm_csv} \
    drep_output

    """

    stub:
    """
    mkdir -p drep_output/dereplicated_genomes
    """
}
