process RENAME_AND_CHECK_SIZE_BINS {

    tag "${bin}"

    input:
    val(name)
    path(bin)

    output:
    path("${name}.${bin.baseName}.fa"), emit: renamed, optional: true

    script:
    """
    SIZE=\$(stat -Lc "%s" ${bin})
    if (( \$SIZE > 50000)) && (( \$SIZE < 20000000)); then
        echo "\${SIZE}"
        cp ${bin} ${name}.${bin.baseName}.fa
    else
        echo "Skipping ${bin} because the bin size \${SIZE} is not between 50kb and 20Mb"
    fi
    """
}
