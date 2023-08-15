process CONSOLIDATE_BINS {

    publishDir(
        path: "${params.outdir}/",
        mode: 'copy',
        failOnError: true
    )

    container 'quay.io/biocontainers/biopython:1.75'

    input:
    path(binner_folders)
    path(stats_files, stageAs: "stats/*")

    output:
    path("consolidated_bins"), emit: consolidated_bins
    path("consolidated_stats.tsv"), emit: consolidated_stats
    path("dereplicated_bins/*"), emit: dereplicated_bins

    script:
    """
    consolidate_bins.py -i ${binner_folders} -s stats -v
    """
}