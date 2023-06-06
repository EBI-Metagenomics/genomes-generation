process SELECT_CONT_SEQUENCES {

    tag "${name} ${bin}"
    container 'quay.io/biocontainers/biopython:1.75'

    input:
    tuple val(name), path(bin), path(summary), path(names)

    output:
    tuple val(name), path("${name}_${bin.baseName}_clean.fa"), emit: cleaned_fasta
    tuple val(name), path("${name}.cont-contigs.txt"), path("*.cont-stats.tsv"), emit: contamination_stats

    script:
    """
    echo "detect contamination"
    detect_contamination.py -s ${summary} -n ${names} -i ${name}

    LEN=\$(wc -l ${name}.cont-contigs.txt)
    if [ \$LEN -eq 0 ]; then
        echo "no contamination"
    else
        echo "cleaning fasta"
        select_seqs_notin_IDs.py -i ${bin} -d ${name}.cont-contigs.txt -o "${name}_${bin.baseName}_clean.fa"
    fi
    """

    stub:
    """
    touch \
        "${name}.cont-contigs.txt" \
        "${name}.cont-stats.tsv" \
        "${name}_${bin.baseName}_clean.fa"
    """
}