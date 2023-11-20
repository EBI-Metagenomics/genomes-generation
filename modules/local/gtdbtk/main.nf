process GTDBTK {

    container 'quay.io/biocontainers/gtdbtk:2.3.0--pyhdfd78af_2'

    containerOptions "${ workflow.containerEngine == 'singularity' ?
        "--bind ${gtdbtk_refdata}:/opt/gtdbtk_refdata":
        "-v ${gtdbtk_refdata}:/opt/gtdbtk_refdata" }"

    input:
    path(genomes_fna, stageAs: "genomes_dir/*")
    path(gtdbtk_refdata)

    output:
    path "gtdbtk_results.tar.gz", emit: gtdbtk_output_tarball
    path "versions.yml"         , emit: versions


    script:
    // TODO: tweak the cpus based on the number of genomes
    """
    export GTDBTK_DATA_PATH=/opt/gtdbtk_refdata

    gtdbtk classify_wf \
    --cpus ${task.cpus} \
    --pplacer_cpus ${task.cpus} \
    --genome_dir genomes_dir \
    --extension fa \
    --skip_ani_screen \
    --out_dir gtdbtk_results

    tar -czf gtdbtk_results.tar.gz gtdbtk_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """

    stub:
    """
    export GTDBTK_DATA_PATH=/opt/gtdbtk_refdata

    mkdir gtdbtk_results

    mkdir -p gtdbtk_results/classify
    touch gtdbtk_results/classify/gtdbtk.bac120.summary.tsv
    touch gtdbtk_results/classify/gtdbtk.ar122.summary.tsv

    tar -czf gtdbtk_results.tar.gz gtdbtk_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
