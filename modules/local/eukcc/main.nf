process LINKTABLE {

    tag "${meta.id} ${binner}"

    // Multi-Package BioContainer
    // TODO: this is using old version of biopython and pysam
    // FIXME: EukCC include biopython and pysam in the EukCC conda package
    container 'quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0'

    input:
    tuple val(meta), path(fasta), path(bam), path(bai), path(bins, stageAs: "bins/*")
    val(binner)

    output:
    tuple val(meta), path("*.links.csv"), emit: links_table
    path "versions.yml"                 , emit: versions

    script:
    """
    mkdir -p bins
    BINS=\$(ls bins | grep -v "unbinned" | wc -l)
    if [ \$BINS -eq 0 ]; then
        echo "creating empty links file"
        touch ${meta.id}.${binner}.links.csv
    else
        binlinks.py --ANI 99 \
        --within 1500 \
        --out ${meta.id}.${binner}.links.csv \
        --debug \
        --bindir bins \
        --bam ${bam[0]}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
        pysam: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """
}

/*
 * EukCC
--- Examples of outputs:
eukcc.csv:
bin     completeness    contamination   ncbi_lng
accession_bin.13.fa        3.88    0.0     1-131567-2759-2698737-33630-2864-89954-252141
accession_merged.0.fa      73.6    0.6     1-131567-2759-33154-4751-451864-5204-452284-1538075-162474-742845

merged.csv:
merged  bins
accession_merged.0.fa      accession_bin.12.fachild:accession_bin.9.fa

merged_bins/
accession_merged.0.fa
*/
process EUKCC {

    tag "${meta.id} ${binner}"

    container 'quay.io/biocontainers/eukcc:2.1.0--pypyhdfd78af_0'

    input:
    val binner
    tuple val(meta), path(links), path(bins, stageAs: "bins/*")
    path eukcc_db

    output:
    tuple val(meta), path("*_merged_bins/${binner}_${meta.id}_merged_bins/*"), optional: true, emit: eukcc_merged_bins
    tuple val(meta), path("${meta.id}_${binner}.eukcc.csv"),                                   emit: eukcc_csv
    tuple val(meta), path("${meta.id}_${binner}.merged_bins.csv"),             optional: true, emit: eukcc_merged_csv
    path "versions.yml",                                                                       emit: versions

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
    echo "EukCC finished"

    cp ${binner}_${meta.id}_merged_bins/eukcc.csv ${meta.id}_${binner}.eukcc.csv
    lines=\$(wc -l < ${binner}_${meta.id}_merged_bins/merged_bins.csv)
    if [ \${lines} -gt 1 ]; then
        cp ${binner}_${meta.id}_merged_bins/merged_bins.csv ${meta.id}_${binner}.merged_bins.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EukCC: \$( eukcc -v | grep -o '[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' )
    END_VERSIONS
    """
}