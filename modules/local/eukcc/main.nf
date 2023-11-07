process LINKTABLE {
    tag "${meta.id} ${binner}"

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    tuple val(meta), path(fasta), path(bam), path(bai), path(bins, stageAs: "bins/*")
    val(binner)

    output:
    tuple val(meta), path("*.links.csv"), emit: links_table

    script:
    """
    BINS=\$(ls bins | grep -v "unbinned" | wc -l)
    if [ \$BINS -eq 0 ]; then
        echo "creating empty links file"
        touch ${meta.id}.${binner}.links.csv
    else
        binlinks.py --ANI 99 \
        --within 1500 \
        --out ${meta.id}.${binner}.links.csv \
        --bindir bins \
        --bam ${bam[0]} -d
    fi
    """
}

/*
 * EukCC
*/
process EUKCC {
    tag "${meta.id} ${binner}"

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    val binner
    tuple val(meta), path(links), path(bins, stageAs: "bins/*")
    path eukcc_db

    output:
    tuple val(meta), path("*_merged_bins"), emit: eukcc_results
    tuple val(meta), path("${meta.id}_${binner}.eukcc.csv"), emit: eukcc_csv
    tuple val(meta), path("${meta.id}_${binner}.merged_bins.csv"), emit: eukcc_merged_csv

    script:
    """
    eukcc folder \
        --improve_percent 10 \
        --n_combine 1 \
        --threads ${task.cpus} \
        --improve_ratio  5 \
        --links ${links} \
        --min_links 100 \
        --suffix .fa \
        --db ${eukcc_db} \
        --out ${binner}_${meta.id}_merged_bins \
        --prefix "${binner}_${meta.id}_merged." \
        bins

    cp *_merged_bins/eukcc.csv ${meta.id}_${binner}.eukcc.csv
    cp *_merged_bins/merged_bins.csv ${meta.id}_${binner}.merged_bins.csv
    """
}