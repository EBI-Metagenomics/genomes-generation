process CHECKM2 {

    tag "${name} ${meta.id}"

    container 'quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0'

    input:
    val(name)
    tuple val(meta), path(bins, stageAs: "bins/*")  // bins can be a list or directory (checkm2 supports both)
    path checkm2_db

    output:
    tuple val(meta), path(bins), path("all_stats.csv"), emit: checkm2_results
    tuple val(meta), path("${name}_filtered_genomes"), optional: true, emit: checkm2_results_filtered
    tuple val(meta), path("${name}_filtered_genomes.tsv"), optional: true, emit: checkm2_results_filtered_stats

    script:
    """
    echo "checkm predict"
    checkm2 predict --threads ${task.cpus} \
        --input bins \
        -x fa \
        --output-directory ${name}_checkm_output \
        --database_path ${checkm2_db}

    echo "checkm table"
    echo "genome,completeness,contamination" > ${name}_checkm2.tsv
    tail -n +2 ${name}_checkm_output/quality_report.tsv | cut -f1-3 | tr '\\t' ',' >> ${name}_checkm2.tsv

    echo "genome,completeness,contamination" > all.stats.csv
    tail -n +2 ${name}_checkm2.tsv | sed 's/\\,/.fa\\,/' >> all_stats.csv

    echo "filter genomes"
    echo "bin\tcompleteness\tcontamination" > ${name}_filtered_genomes.tsv
    cat ${name}_checkm2.tsv | \
        tr ',' '\\t' |\
        grep -v "completeness" |\
        awk '{{if(\$2>=50 && \$2<=100 && \$3>=0 && \$3<=5){{print \$0}}}}' >> ${name}_filtered_genomes.tsv

    echo "choose genomes"
    mkdir -p ${name}_filtered_genomes
    for i in \$(cat ${name}_filtered_genomes.tsv | grep -v "completeness" | cut -f1 ); do
        cp bins/\${i}.* ${name}_filtered_genomes
    done
    """
}