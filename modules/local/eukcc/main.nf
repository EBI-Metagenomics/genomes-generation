/*
 * EukCC
--- Examples of outputs:
eukcc.tsv (published as <bins>.eukcc.csv, tab separated):
bin     completeness    contamination   ncbi_lng        ncbi_lng_names
accession_bin.13.fa        3.88    0.0     1-131567-2759-2698737-33630-2864-89954-252141  ...
accession_merged.0.fa      73.6    0.6     1-131567-2759-33154-4751-451864-5204-452284-1538075-162474-742845  ...
accession_bin.7.fa         NA      NA      NA      NA   <- bins that could not be estimated (no proteins / no marker set)

merged_bins.tsv (published as <bins>.merged_bins.csv):
merged  bins
accession_merged.0.fa      accession_bin.12.fachild:accession_bin.9.fa

merged_bins/
accession_merged.0.fa
*/
process EUKCC {

    label 'process_long'
    tag "${meta.id} ${bins.baseName}"

    container 'quay.io/biocontainers/eukcc:2.3.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(links), path(bins)
    path eukcc_db

    output:
    tuple val(meta), path("${bins.baseName}_results/merged_bins"),        emit: eukcc_merged_bins
    tuple val(meta), path("${bins.baseName}.eukcc.csv"),                  emit: eukcc_csv
    tuple val(meta), path("${bins.baseName}.merged_bins.csv"),            emit: eukcc_merged_csv
    path "versions.yml",                                                  emit: versions
    path "progress.log",                                                  emit: progress_log

    script:
    """
    mkdir -p ${bins.baseName}_results/merged_bins

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EukCC: \$( eukcc -v | grep -o '[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' )
    END_VERSIONS

    export BINS=\$(ls $bins | wc -l)
    if [ \$BINS -eq 0 ]; then
        echo "No bins in input"
        touch ${bins.baseName}.eukcc.csv ${bins.baseName}.merged_bins.csv
    else
        eukcc folder \
            --improve_percent 10 \
            --n_combine 1 \
            --threads ${task.cpus} \
            --improve_ratio  5 \
            --links ${links} \
            --min_links 100 \
            --db ${eukcc_db} \
            --out ${bins.baseName}_results \
            --prefix "${bins.baseName}_merged." \
            $bins

        mv ${bins.baseName}_results/eukcc.tsv ${bins.baseName}.eukcc.csv
        mv ${bins.baseName}_results/merged_bins.tsv ${bins.baseName}.merged_bins.csv
    fi

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}\t${bins.baseName}
        bins: \$BINS, merged: \$(ls ${bins.baseName}_results/merged_bins | wc -l)
    END_LOGGING
    """
}
