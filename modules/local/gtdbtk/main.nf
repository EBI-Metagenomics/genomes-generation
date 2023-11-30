process GTDBTK {

    container 'quay.io/microbiome-informatics/gtdb-tk:2.3.0'

    containerOptions "${ workflow.containerEngine == 'singularity' ?
        "--bind ${gtdbtk_refdata}:/opt/gtdbtk_refdata":
        "-v ${gtdbtk_refdata}:/opt/gtdbtk_refdata" }"

    input:
    path(genomes_fna, stageAs: "genomes_dir/*")
    path(gtdbtk_refdata)

    output:
    path "gtdbtk_results.tar.gz", emit: gtdbtk_output_tarball
    path "ncbi_taxonomy.txt", emit: ncbi_taxonomy
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

    echo "Create NCBI taxonomy"
    gtdb_to_ncbi_majority_vote.py \
        --gtdbtk_output_dir gtdbtk_results \
        --output_file ncbi_taxonomy.txt \
        --ar53_metadata_file ${gtdbtk_refdata}/ar53_metadata_r214.tsv \
        --bac120_metadata_file ${gtdbtk_refdata}/bac120_metadata_r214.tsv

    echo "Compress GTDB-Tk"
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
    touch ncbi_taxonomy.txt

    tar -czf gtdbtk_results.tar.gz gtdbtk_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
