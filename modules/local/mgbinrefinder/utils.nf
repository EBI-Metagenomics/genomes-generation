process RENAME_AND_CHECK_SIZE_BINS {

    label 'process_low'
    tag "${meta.id}"

    input:
    val(name)
    tuple val(meta), path(bins, stageAs: "bins_dir/*")

    output:
    tuple val(meta), path("out/*"), emit: renamed, optional: true
    path "progress.log"           , emit: progress_log

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
    cd ..

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}\t${name}
        bins_dir: \$(ls bins_dir | wc -l), filtered_out: \$(ls out | wc -l)
    END_LOGGING
    """
}
