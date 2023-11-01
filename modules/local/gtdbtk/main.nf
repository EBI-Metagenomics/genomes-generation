process GTDBTK {

    container 'quay.io/microbiome-informatics/gtdb-tk:2.1.0'

    containerOptions "${ workflow.containerEngine == 'singularity' ?
        "--bind ${gtdbtk_refdata}:/opt/gtdbtk_refdata":
        "-v ${gtdbtk_refdata}:/opt/gtdbtk_refdata" }"

    input:
    path(genomes_fna, stageAs: "genomes_dir/*")
    path(gtdbtk_refdata)

    output:
    path 'gtdbtk_results.tar.gz', emit: gtdbtk_output_tarball

    script:
    // TODO: tweak the cpus based on the number of genomes
    """
    GTDBTK_DATA_PATH=/opt/gtdbtk_refdata \
    gtdbtk classify_wf \
    --cpus ${task.cpus} \
    --pplacer_cpus ${task.cpus} \
    --genome_dir genomes_dir \
    --extension fa \
    --out_dir gtdbtk_results

    tar -czf gtdbtk_results.tar.gz gtdbtk_results
    """

    stub:
    """
    mkdir gtdbtk_results

    mkdir -p gtdbtk_results/classify
    touch gtdbtk_results/classify/gtdbtk.bac120.summary.tsv
    touch gtdbtk_results/classify/gtdbtk.ar122.summary.tsv

    tar -czf gtdbtk_results.tar.gz gtdbtk_results
    """
}
