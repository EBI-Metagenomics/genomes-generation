process LINKTABLE {
    tag "${name} ${binner}"
    publishDir(
        path: "${params.outdir}/eukcc/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/microbiome-informatics/eukrecover.python3base:v1'

    input:
    tuple val(name), path(bam), path(bindir)
    val(binner)

    output:
    tuple val(name), path("*.links.csv"), path(bindir), emit: links_table

    script:
    """
    BINS=\$(ls ${bindir} | grep -v "unbinned" | wc -l)
    if [ \$BINS -eq 0 ]; then
        echo "creating empty links file"
        touch ${name}.${binner}.links.csv
    else
        binlinks.py  --ANI 99 --within 1500 --out ${name}.${binner}.links.csv --bindir ${bindir} --bam ${bam[0]}
    fi
    """
}

/*
 * EukCC
*/
process EUKCC {
    tag "${name} ${binner}"

    publishDir(
        path: "${params.outdir}/eukcc/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    val binner
    tuple val(name), path(links), path(bindir)
    path eukcc_db

    output:
    tuple val(name), path("*_merged_bins"), emit: eukcc_results
    tuple val(name), path("${name}_${binner}.eukcc.csv"), emit: eukcc_csv

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
        --out ${binner}_${name}_merged_bins \
        --prefix "${binner}_${name}_merged." \
        ${bindir}

    cp *_merged_bins/eukcc.csv ${name}_${binner}.eukcc.csv
    """
}

process EUKCC_MAG {
    tag "${name}"

    publishDir(
        path: "${params.outdir}/eukcc_mags/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    tuple val(name), path(genomes_list, stageAs: "mags_dir/*")
    path eukcc_db

    output:
    tuple val(name), path("output/*eukcc.csv"), emit: eukcc_mags_results

    script:
    """
    eukcc --debug folder \
            --threads ${task.cpus} \
            --db ${eukcc_db} \
            --out output \
            mags_dir
    """

    stub:
    """
    touch mags.eukcc.csv
    """
}

process EUKCC_SINGLE {
    tag "${genome}"

    publishDir(
        path: "${params.outdir}/eukcc_mags/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/microbiome-informatics/eukcc:latest'

    input:
    path(genome)
    path eukcc_db

    output:
    path("*eukcc.csv"), emit: eukcc_mags_results

    script:
    """
    eukcc --debug single \
            --threads ${task.cpus} \
            --db ${eukcc_db} \
            --out output \
            ${genome}
    cp output/eukcc.csv ${genome.baseName}.eukcc.csv
    """

    stub:
    """
    touch mags.eukcc.csv
    """
}