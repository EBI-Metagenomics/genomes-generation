process CAT {
    tag "${meta.id} ${bin_fa}"

    container 'quay.io/microbiome-informatics/cat:5.2.3'

    input:
    tuple val(meta), path(bin_fa)
    path cat_db
    path taxonomy_db
    path diamond_db

    output:
    tuple val(meta), path(bin_fa), path('*.summary.txt'), path('*.contig2classification.official_names.txt'), emit: cat_results

    script:
    """
    echo "[MAG clean-up] Analysing contigs"
    CAT contigs -n ${task.cpus} -c ${bin_fa} \
    --path_to_diamond ${diamond_db} -d ${cat_db} -t ${taxonomy_db} --out_prefix ${meta.id}

    echo "[MAG clean-up] Adding taxonomy names"
    CAT add_names \
    -i ${meta.id}.contig2classification.txt -o ${meta.id}.contig2classification.official_names.txt \
    -t ${taxonomy_db} --only_official

    echo "[MAG clean-up] Summarizing output"
    CAT summarise -c ${bin_fa} \
    -i ${meta.id}.contig2classification.official_names.txt -o ${meta.id}.summary.txt
    """

    stub:
    """
    touch ${meta.id}.contig2classification.official_names.txt ${meta.id}.summary.txt
    """
}