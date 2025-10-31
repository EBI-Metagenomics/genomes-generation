process CHECKM2 {

    label 'process_medium'

    tag "${meta.id}"

    container 'quay.io/biocontainers/checkm2:1.0.1--pyh7cba7a3_0'

    input:
    tuple val(meta), path(bins, stageAs: "bins/*")
    path checkm2_db

    output:
    tuple val(meta), path("bins_folder"), path("${prefix}_all_stats.csv")   , emit: bins_and_stats
    // TODO filtered_genomes and filtered_stats are not used in the pipeline
    // consider removing them or replacing the module with nf-core/checkm2
    tuple val(meta), path("${prefix}_filtered_genomes")                     , emit: filtered_genomes
    tuple val(meta), path("${prefix}_filtered_genomes.tsv")                 , emit: filtered_stats
    path "versions.yml"                                                     , emit: versions
    path "progress.log"                                                     , emit: progress_log

    // NOTE:
    // Checkm2 works with list of files OR folder as --input
    // bins can be a folder with bins OR folder with folder of bins
    // there is a check for both structures
    // bins finally moved to bins_folder

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS

    mkdir -p bins_folder "${prefix}_filtered_genomes"
    touch "${prefix}_all_stats.csv" "${prefix}_filtered_genomes.tsv"

    set +e

    export BINS=\$(restructure_input.py -i bins -o bins_folder)

    echo "Check the number of bins"
    if [ \$BINS -eq 0 ]; then
        echo "Bins folder is empty"
        exit 0
    fi

    echo "checkm predict"
    mkdir -p checkm_tmp
    checkm2 predict --threads ${task.cpus} \
        --input bins_folder \
        -x fa \
        --output-directory ${prefix}_checkm_output \
        --database_path ${checkm2_db} \
        --tmpdir checkm_tmp
    CHECMK2_EXITCODE="\$?"

    if [ "\$CHECMK2_EXITCODE" == "1" ]; then
        echo "Checkm2 exit code \$CHECMK2_EXITCODE"
        if [ ! -e "${prefix}_checkm_output/checkm2.log" ]; then
            echo "checkm2.log does not exist. Exit"
            echo "Error" >&2
            exit \$CHECMK2_EXITCODE
        else
            echo "checkm2.log exists -> checking if DIAMOND failed"
            if grep -q "No DIAMOND annotation was generated. Exiting" "${prefix}_checkm_output/checkm2.log"; then
                echo "No DIAMOND annotation was generated"
                touch "${prefix}_all_stats.csv" "${prefix}_filtered_genomes.tsv"
                mkdir "${prefix}_filtered_genomes"
                exit 0
            else
                echo "It is not DIAMOND, sorry, check manually. Exit" >&2
                exit \$CHECMK2_EXITCODE
            fi
        fi
    fi
    set -e

    echo "checkm table"
    echo "genome,completeness,contamination" > ${prefix}_checkm2.tsv
    tail -n +2 ${prefix}_checkm_output/quality_report.tsv | cut -f1-3 | tr '\\t' ',' >> ${prefix}_checkm2.tsv

    awk -F, 'NR == 1 {print; next} {OFS=","; \$1 = \$1 ".fa"; print}' ${prefix}_checkm2.tsv > ${prefix}_all_stats.csv

    echo "filter genomes"
    echo "bin\tcompleteness\tcontamination" > ${prefix}_filtered_genomes.tsv
    cat ${prefix}_checkm2.tsv | \
        tr ',' '\\t' |\
        grep -v "completeness" |\
        awk '{{if(\$2>=50 && \$2<=100 && \$3>=0 && \$3<=5){{print \$0}}}}' >> ${prefix}_filtered_genomes.tsv

    echo "choose genomes"
    mkdir -p ${prefix}_filtered_genomes
    for i in \$(cat ${prefix}_filtered_genomes.tsv | grep -v "completeness" | cut -f1 ); do
        cp bins_folder/\${i}.* ${prefix}_filtered_genomes
    done

    cat <<-END_LOGGING >> progress.log
    ${meta.id}\t${task.process}
        bins_folder: \$(ls bins_folder | wc -l), filtered: \$(ls ${prefix}_filtered_genomes | wc -l)
    END_LOGGING
    """
}
