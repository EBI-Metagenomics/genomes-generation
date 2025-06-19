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

--- Exit codes workaround:
if exit code is 204 - metaeuk return empty faa -> no results for eukcc
*/
process EUKCC {

    label 'process_long'
    tag "${meta.id} ${bins.baseName}"

    container 'quay.io/microbiome-informatics/eukcc:2.1.3'

    input:
    tuple val(meta), path(links), path(bins)
    path eukcc_db

    output:
    tuple val(meta), path("${meta.id}_merged_bins/merged_bins"),    emit: eukcc_merged_bins
    tuple val(meta), path("${meta.id}.eukcc.csv"),                  emit: eukcc_csv
    tuple val(meta), path("${meta.id}.merged_bins.csv"),            emit: eukcc_merged_csv
    path "versions.yml",                                                      emit: versions
    path "progress.log",                                                      emit: progress_log

    script:
    """
    mkdir -p bins ${meta.id}_merged_bins/merged_bins
    touch ${meta.id}.eukcc.csv ${meta.id}.merged_bins.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EukCC: \$( eukcc -v | grep -o '[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' )
    END_VERSIONS

    export BINS=\$(ls bins | wc -l)
    if [ \$BINS -eq 0 ]; then
        echo "No bins in input"
        mkdir -p ${meta.id}_merged_bins/merged_bins
        cat <<-END_LOGGING > progress.log
        ${meta.id}\t${task.process}\t${bins.baseName}
            bins: 0, merged: 0
    END_LOGGING
    else
        set +e

        eukcc folder \
            --improve_percent 10 \
            --n_combine 1 \
            --threads ${task.cpus} \
            --improve_ratio  5 \
            --links ${links} \
            --min_links 100 \
            --suffix .fa \
            --db ${eukcc_db} \
            --out ${meta.id}_merged_bins \
            --prefix "${meta.id}_merged." \
            bins
        EUKCC_EXITCODE="\$?"

        if [ "\$EUKCC_EXITCODE" == "0" ]; then
            echo "EukCC finished"
            mv ${meta.id}_merged_bins/eukcc.csv ${meta.id}.eukcc.csv
            mv ${meta.id}_merged_bins/merged_bins.csv ${meta.id}.merged_bins.csv
        fi

        if [ "\$EUKCC_EXITCODE" == "204" ]; then
            echo "Metaeuk returned zero proteins"
        fi

        mkdir -p ${meta.id}_merged_bins/merged_bins

        set -e

        if [ "\$EUKCC_EXITCODE" == "0" ] || [ "\$EUKCC_EXITCODE" == "204" ]; then
            cat <<-END_LOGGING > progress.log
            ${meta.id}\t${task.process}\t${bins.baseName}
                bins: \$(ls bins | wc -l), merged: \$(ls ${meta.id}_merged_bins/merged_bins | wc -l)
    END_LOGGING
            exit 0
        else:
            exit \$EUKCC_EXITCODE
        fi

    fi
    """
}
