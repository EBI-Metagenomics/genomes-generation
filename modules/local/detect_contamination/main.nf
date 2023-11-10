process DETECT_CONTAMINATION {

    tag "${meta.id} ${bin}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(bin), path(summary), path(names)

    output:
    tuple val(meta), path("${meta.id}_${bin.baseName}_clean.fa"),                                                       emit: cleaned_fasta
    tuple val(meta), path("${meta.id}.contamination_contigs.txt"), path("*.contamination_contigs.tsv"), optional: true, emit: contamination_stats
    path "versions.yml"                                                                                               , emit: versions

    // TODO: merge detect_decontamination and select_seqs_notin_ids into a single python script
    script:
    """
    echo "detect contamination"
    detect_contamination.py -s ${summary} -n ${names} -i ${meta.id}

    lines_count=\$(wc -l < "${meta.id}.contamination_contigs.txt")
    if [ "\$lines_count" -eq 0 ]; then
        echo "no contamination"
        cp ${bin} ${meta.id}_${bin.baseName}_clean.fa
    else
        echo "cleaning fasta"
        select_seqs_notin_ids.py -i ${bin} -d ${meta.id}.contamination_contigs.txt -o "${meta.id}_${bin.baseName}_clean.fa"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('biopython').version)")
    END_VERSIONS
    """

    stub:
    """
    touch \
        "${meta.id}.contamination_contigs.txt" \
        "${meta.id}.cont-stats.tsv" \
        "${meta.id}_${bin.baseName}_clean.fa"
    """
}