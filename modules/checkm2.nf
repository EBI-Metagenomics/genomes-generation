process CHECKM2 {

    publishDir(
        path: "${params.outdir}/checkm2",
        mode: 'copy',
        failOnError: true
    )

    input:
    path bins
    path checkm_db

    output:
    path "all.stats.clean", emit: checkm_output

    script:
    """
    echo "checkm predict"
    checkm2 predict --threads ${task.cpus} --input ${bins} -x fa --output-directory checkm_output

    echo "checkm table"
    echo "genome,completeness,contamination" > all.stats.clean
    tail -n +2 checkm_output/quality_report.tsv | cut -f1-3 | tr '\t' ',' | sed 's/\\,/\\.fa\\,/' >> all.stats.clean
    """
}