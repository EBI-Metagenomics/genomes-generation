process BUSCO {
    label 'process_medium'
    tag "${bin}"

    container 'quay.io/biocontainers/busco:5.4.7--pyhdfd78af_0'

    input:
    path bin
    path busco_db

    output:
    path "short_summary.specific*.txt", emit: busco_summary
    path "versions.yml"                , emit: versions

    script:
    """
    busco  --offline \
            -i ${bin} \
            -m 'genome' \
            -o out \
            --auto-lineage-euk \
            --download_path ${busco_db} \
            -c ${task.cpus}

    mv out/short_summary.specific*.out.txt "short_summary.specific_${bin.baseName}.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
