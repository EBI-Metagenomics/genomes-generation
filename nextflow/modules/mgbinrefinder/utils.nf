process RENAME_AND_CHECK_SIZE_BINS {

    tag "${bin}"

    input:
    val(name)
    tuple val(meta), path(bin)

    output:
    tuple val(meta), path("out/*"), emit: renamed, optional: true

    script:
    """
    SIZE=\$(stat -Lc "%s" ${bin})
    mkdir -p out
    if (( \$SIZE > 50000)) && (( \$SIZE < 20000000)); then
        echo "\${SIZE}"
        cp ${bin} out
    else
        echo "Skipping ${bin} because the bin size \${SIZE} is not between 50kb and 20Mb"
    fi
    """
}