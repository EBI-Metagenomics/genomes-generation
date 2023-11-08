process RENAME_AND_CHECK_SIZE_BINS {

    tag "${meta.id}"

    input:
    val(name)
    tuple val(meta), path(bins, stageAs: "bins_dir/*")

    output:
    tuple val(meta), path("out/*"), emit: renamed, optional: true

    script:
    """
    mkdir -p out
    for bin in \$(find bins_dir -type f); do
        SIZE=\$(stat -c "%s" \${bin})
        if (( \$SIZE > 50000)) && (( \$SIZE < 20000000)); then
            echo "\${SIZE}"
            cp \${bin} out
        else
            echo "Skipping \${bin} because the bin size \${SIZE} is not between 50kb and 20Mb"
        fi
    done
    """
}
