process CAT {
    label 'process_medium'
    tag "${meta.id} ${bin_fa}"

    container 'quay.io/microbiome-informatics/cat:5.2.3'

    input:
    tuple val(meta), path(bin_fa)
    path cat_db_folder
    path taxonomy_db
    path cat_diamond_db

    output:
    tuple val(meta), path(bin_fa), path('*.summary.txt'), path('*.contig2classification.official_names.txt'), emit: cat_results
    path "versions.yml"                                                                                     , emit: versions

    script:
    """
    echo "[MAG clean-up] Analysing contigs"
    mkdir -p cat_tmp
    CAT contigs \
    -n ${task.cpus} \
    -c ${bin_fa} \
    -d ${cat_db_folder} \
    -t ${taxonomy_db} \
    --path_to_diamond ${cat_diamond_db} \
    --tmpdir cat_tmp \
    --out_prefix ${meta.id}

    echo "[MAG clean-up] Adding taxonomy names"
    CAT add_names \
    -i ${meta.id}.contig2classification.txt -o ${meta.id}.contig2classification.official_names.txt \
    -t ${taxonomy_db} --only_official

    echo "[MAG clean-up] Summarizing output"
    CAT summarise -c ${bin_fa} \
    -i ${meta.id}.contig2classification.official_names.txt -o ${meta.id}.summary.txt

    echo "Remove diamond alignment"
    rm -f "${meta.id}.alignment.diamond"

    echo "Remove predicted proteins"
    rm -f "${meta.id}.predicted_proteins.faa"
    rm -f "${meta.id}.predicted_proteins.gff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CAT: \$(CAT --version | sed "s/CAT v//; s/(.*//")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.contig2classification.official_names.txt ${meta.id}.summary.txt
    """
}