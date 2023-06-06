process CAT {
    tag "${name} ${bin_fa}"

    container 'quay.io/microbiome-informatics/cat:5.2.3'

    publishDir(
        path: "${params.outdir}/cleaning_bins/",
        mode: 'copy',
        failOnError: true
    )

    input:
    tuple val(name), path(bin_fa)
    path cat_db
    path taxonomy_db
    path diamond_db

    output:
    tuple val(name), path(bin_fa), path('*.summary.txt'), path('*.contig2classification.official_names.txt'), emit: cat_results

    script:
    """
    echo "[MAG clean-up] Analysing contigs"
    CAT contigs -n ${task.cpus} -c ${bin_fa} \
    --path_to_diamond ${diamond_db} -d ${cat_db} -t ${taxonomy_db} --out_prefix ${name}

    echo "[MAG clean-up] Adding taxonomy names"
    CAT add_names \
    -i ${name}.contig2classification.txt -o ${name}.contig2classification.official_names.txt \
    -t ${taxonomy_db} --only_official

    echo "[MAG clean-up] Summarizing output"
    CAT summarise -c ${bin_fa} \
    -i ${name}.contig2classification.official_names.txt -o ${name}.summary.txt
    """

    stub:
    """
    touch ${name}.contig2classification.official_names.txt ${name}.summary.txt
    """
}