process GUNC {
    label 'process_low'
    tag "${meta.id} ${fasta}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.gunc:v4'

    input:
    tuple val(meta), path(fasta)
    file gunc_db

    output:
    tuple val(meta), path("${fasta.baseName.replaceAll('_clean', '')}.${fasta.extension}"), path('*_gunc_*'), emit: tuple_gunc_result
    path 'gunc_contaminated.txt'                                                                            , emit: gunc_result
    path "versions.yml"                                                                                     , emit: versions

    script:
    """
    gunc run -t ${task.cpus} \
    -i ${fasta} \
    --db_file ${gunc_db}

    ### gunc contaminated genomes ###
    cat GUNC.*.maxCSS_level.tsv | grep -v "n_genes_mapped" | awk '{if(\$8 > 0.45 && \$9 > 0.05 && \$12 > 0.5) print\$1}' > gunc_contaminated.txt

    # gunc_contaminated.txt could be empty - that means genome is OK
    # gunc_contaminated.txt could have this genome inside - that means gunc filtered this genome

    if [ -s gunc_contaminated.txt ]; then
        touch ${fasta.baseName}_gunc_empty.txt
    else
        touch ${fasta.baseName}_gunc_complete.txt
    fi

    mv "$fasta" "${fasta.baseName.replaceAll('_clean', '')}.${fasta.extension}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunc: \$( gunc --version )
    END_VERSIONS
    """

    stub:
    """
    touch gunc_contaminated.txt
    touch ${fasta.baseName}_gunc_empty.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunc: \$( gunc --version )
    END_VERSIONS
    """
}
