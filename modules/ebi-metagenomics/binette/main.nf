process BINETTE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/binette:1.1.2--pyh7e72e81_0'
        : 'biocontainers/binette:1.1.2--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(input_binning, stageAs: "input_binning/*"), path(contigs), path(proteins)
    val input_type
    path checkm2_db

    output:
    tuple val(meta), path("${meta.id}_final_bins/*")                  , emit: refined_bins
    tuple val(meta), path("${meta.id}_final_bins_quality_reports.tsv"), emit: refined_bins_report
    path "progress.log"                                               , emit: progress_log
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def protein_arg = proteins ? "--proteins ${proteins}" : ""
    def input_arg   = ""

    if (input_type == 'fasta') {
        input_arg = "--bin_dirs"
    } else if (input_type == 'tsv') {
        input_arg = "--contig2bin_tables"
    } else {
        error "Invalid input_type: ${input_type}. Must be 'fasta' or 'tsv'"
    }

    """
    binette \\
        ${input_arg} input_binning/* \\
        ${protein_arg} \\
        --contigs ${contigs} \\
        --threads ${task.cpus} \\
        --checkm2_db ${checkm2_db} \\
        -o . \\
        ${args}

    # add run accession to all bin names
    for file in final_bins/*.fa*; do mv "\$file" "final_bins/${meta.id}_\$(basename "\$file")"; done
    # add run accession to the bin folder and report with stats
    mv final_bins "${meta.id}_final_bins"
    mv final_bins_quality_reports.tsv ${meta.id}_final_bins_quality_reports.tsv

    # Count number of bins before and after refinement
    BINNER_COUNTS=""
    BINNER_INDEX=1
    for binner_dir in input_binning/*; do
        if [ -d "\$binner_dir" ]; then
            count=\$(ls -1 "\$binner_dir" | wc -l)
            BINNER_COUNTS="\${BINNER_COUNTS}binner\${BINNER_INDEX}: \${count}, "
            BINNER_INDEX=\$((BINNER_INDEX + 1))
        fi
    done

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}
        \${BINNER_COUNTS}refined: \$(ls -1 ${meta.id}_final_bins/ | wc -l)
    END_LOGGING

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$( binette --version )
    END_VERSIONS
    """

    stub:
    """
    mkdir ${meta.id}_final_bins
    touch ${meta.id}_final_bins/bin_1.fasta
    touch ${meta.id}_final_bins_quality_reports.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$( binette --version )
    END_VERSIONS
    """
}
