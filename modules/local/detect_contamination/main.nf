process DETECT_CONTAMINATION {

    tag "${meta.id} ${bin}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(bin), path(summary), path(names)

    output:
    // TODO: rename .cont-contigs. -> _contamination_contigs.
    tuple val(meta), path("${meta.id}_${bin.baseName}_clean.fa"), emit: cleaned_fasta
    tuple val(meta), path("${meta.id}.cont-contigs.txt"), path("*.cont-contigs.tsv"), optional: true, emit: contamination_stats

    // TODO: merge detect_decontamination and select_seqs_notin_ids into a single python script
    script:
    """
    echo "detect contamination"
    detect_contamination.py -s ${summary} -n ${names} -i ${meta.id}

    lines_count=\$(wc -l < "${meta.id}.cont-contigs.txt")
    if [ "\$lines_count" -eq 0 ]; then
        echo "no contamination"
        cp ${bin} ${meta.id}_${bin.baseName}_clean.fa
    else
        echo "cleaning fasta"
        select_seqs_notin_ids.py -i ${bin} -d ${meta.id}.cont-contigs.txt -o "${meta.id}_${bin.baseName}_clean.fa"
    fi
    """

    stub:
    """
    touch \
        "${meta.id}.cont-contigs.txt" \
        "${meta.id}.cont-stats.tsv" \
        "${meta.id}_${bin.baseName}_clean.fa"
    """
}