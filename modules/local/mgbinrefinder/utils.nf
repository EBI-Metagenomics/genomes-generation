process RENAME_AND_CHECK_SIZE_BINS {

    label 'process_low'
    tag "${meta.id}"

    input:
    val(name)
    tuple val(meta), path(bins, stageAs: "bins_dir/*")

    output:
    tuple val(meta), path("out/*"), emit: renamed, optional: true

    script:
    """
    mkdir -p out bins_dir
    cd bins_dir
    for bin in \$(ls . ); do
        SIZE=\$(stat -L -c "%s" \${bin})
        if (( \$SIZE > 50000)) && (( \$SIZE < 20000000)); then
            echo "\${SIZE}"
            cp \${bin} ../out
        else
            echo "Skipping \${bin} because the bin size \${SIZE} is not between 50kb and 20Mb"
        fi
    done
    """
}
