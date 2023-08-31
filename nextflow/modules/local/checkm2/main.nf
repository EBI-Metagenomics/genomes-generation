process CHECKM2 {

    tag "${name} ${meta.id}"

    container 'quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0'

    input:
    val(name)
    tuple val(meta), path(bins)
    path checkm_db

    output:
    tuple val(meta), path(bins), path("all.stats.csv"), emit: checkm2_results
    tuple val(meta), path("${meta.id}_${name}_filtered_genomes"), optional: true, emit: checkm2_results_filtered
    tuple val(meta), path("${meta.id}_${name}_filtered_genomes.tsv"), optional: true, emit: checkm2_results_filtered_stats

    script:
    """
    echo "checkm predict"
    checkm2 predict --threads ${task.cpus} \
                    --input ${bins} \
                    -x fa \
                    --output-directory ${name}_checkm_output \
                    --database_path ${checkm_db}

    echo "checkm table"
    echo "genome,completeness,contamination" > ${meta.id}.${name}.checkm2.tsv
    tail -n +2 ${name}_checkm_output/quality_report.tsv | cut -f1-3 | tr '\\t' ',' >> ${meta.id}.${name}.checkm2.tsv

    echo "genome,completeness,contamination" > all.stats.csv
    tail -n +2 ${meta.id}.${name}.checkm2.tsv | sed 's/\\,/.fa\\,/'  >> all.stats.csv

    echo "filter genomes"
    echo "bin\tcompleteness\tcontamination" > ${meta.id}_${name}_filtered_genomes.tsv
    cat ${meta.id}.${name}.checkm2.tsv | \
        tr ',' '\\t' |\
        grep -v "completeness" |\
        awk '{{if(\$2>=50 && \$2<=100 && \$3>=0 && \$3<=5){{print \$0}}}}' >> ${meta.id}_${name}_filtered_genomes.tsv

    echo "choose genomes"
    mkdir -p ${meta.id}_${name}_filtered_genomes
    for i in \$(cat ${meta.id}_${name}_filtered_genomes.tsv | grep -v "completeness" | cut -f1 ); do
        cp \${i}* ${meta.id}_${name}_filtered_genomes
    done
    """
}