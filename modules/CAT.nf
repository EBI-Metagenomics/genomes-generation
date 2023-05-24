process MAG_CLEANUP_CAT {

    container 'quay.io/microbiome-informatics/cat:5.2.3'

    publishDir(
        path: "${params.outdir}/cleaning_bins/",
        mode: 'copy',
        failOnError: true
    )

    input:
    tuple val(sample_name), path(bin_fa)
    path cat_db
    path taxonomy_db
    path diamond_db

    output:
    tuple val(sample_name), path('*.summary.txt'), emit: cat_summary
    tuple val(sample_name), path('*.contig2classification.official_names.txt'), emit: cat_names

    script:
    """
    echo "[MAG clean-up] Analysing contigs"
    CAT contigs -n ${task.cpus} -c ${bin_fa} \
    --path_to_diamond ${diamond_db} -d ${cat_db} -t ${taxonomy_db} --out_prefix ${sample_name}

    echo "[MAG clean-up] Adding taxonomy names"
    CAT add_names \
    -i ${sample_name}.contig2classification.txt -o ${sample_name}.contig2classification.official_names.txt \
    -t ${taxonomy_db} --only_official

    echo "[MAG clean-up] Summarizing output"
    CAT summarise -c ${bin_fa} \
    -i ${sample_name}.contig2classification.official_names.txt -o ${sample_name}.summary.txt
    """
}