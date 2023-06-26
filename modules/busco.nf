process BUSCO {

    publishDir(
        path: "${params.outdir}/busco/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/biocontainers/busco:5.4.7--pyhdfd78af_0'

    input:
    path bin
    path busco_db

    output:
    path "*.short_summary.specific.txt", emit: busco_summary

    script:
    """
    busco  --offline \
              -i ${bin} \
              -m 'genome' \
              -o out \
              --auto-lineage-euk \
              --download_path ${busco_db} \
              -c {task.cpus}

    cp out/short_summary.specific*.out.txt "${bin.baseName}.short_summary.specific.txt"
    """
}