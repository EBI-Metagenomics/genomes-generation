/*
    ~~~~~~~~~~~~~~~~~~
     uncompress
    ~~~~~~~~~~~~~~~~~~
*/
process GUNZIP {

    label 'process_low'
    tag "$meta.id"

    input:
    tuple val(meta), path(compressed_file)

    output:
    tuple val(meta), path("out/*"), emit: uncompressed

    script:
    """
    mkdir out
    cp ${compressed_file} "out/${meta.id}.fasta.gz"
    cd out
    gunzip *
    """
}


process FINALIZE_LOGGING {

    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:0.24.1':
        'quay.io/biocontainers/pandas:0.24.1' }"

    input:
    path(logging_file)
    val(output)

    output:
    path("${output}"), emit: structured_logging

    script:
    """
    logging_stats.py -i ${logging_file} -o ${output}
    """
}


process FILTER_QUALITY {
    tag "${meta.id}"

    label 'process_light'

    input:
    tuple val(meta), path(quality_file), path(bins, stageAs: "input_bins/")
    val delimiter

    output:
    tuple val(meta), path("output_genomes/*"), path("quality_file.csv"), emit: qs50_filtered_genomes, optional: true
    path "progress.log"                                                , emit: progress_log

    script:
    """
    mkdir -p output_genomes
    touch quality_file.csv

    echo "Prepare drep quality"
    grep -v "completeness" ${quality_file} |\
    awk -F "${delimiter}" '{{if(\$2>=50 && \$2<=100 && \$3>=0 && \$3<=5){{print \$0}}}}' |\
    sort -k 2,3 -n | cut -f1 > filtered_genomes.txt || true

    echo "bins count"
    export BINS=\$(cat filtered_genomes.txt | wc -l)
    echo "\$BINS"
    if [ \$BINS -lt 2 ];
    then
        echo "No genomes"
    else
        for i in \$(ls input_bins | grep -w -f filtered_genomes.txt); do
            mv input_bins/\${i} output_genomes; done

        echo "genome,completeness,contamination" > quality_file.csv
        grep -w -f filtered_genomes.txt ${quality_file} | cut -f1-3 | tr '\\t' ',' >> quality_file.csv
    fi

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}
        input_bins: \$(ls input_bins | wc -l)
    END_LOGGING
    """
}
