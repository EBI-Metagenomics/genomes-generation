/*
 * dRep, this workflow dereplicates a set of genomes.
*/
process DREP {
    label 'process_medium'
    tag "${meta.id} ${drep_params}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/drep:3.2.2--pyhdfd78af_0':
        'quay.io/biocontainers/drep:3.2.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(genomes_list, stageAs: "genomes_dir/*"), path(quality_csv)
    val drep_params

    output:
    tuple val(meta), path("drep_output/dereplicated_genomes/*"), optional: true, emit: dereplicated_genomes
    tuple val(meta), path("dereplicated_genomes.txt")          , optional: true, emit: dereplicated_genomes_list
    path "versions.yml"                                                        , emit: versions
    path "progress.log"                                                        , emit: progress_log

    // NOTE:
    // genomes_dir can be a folder with genomes OR folder with folder of genomes
    // there is a check for both structures
    // genomes finally moved to genomes_folder

    script:
    """
    mkdir -p drep_output/dereplicated_genomes genomes_folder
    touch dereplicated_genomes.txt

    echo "Moving genomes to genomes_folder"
    export LEN=\$(restructure_input.py -i genomes_dir -o genomes_folder)

    echo "Checking number of given genomes"
    if [ \$LEN -eq 0 ]; then
        echo "no genomes"
    elif [ \$LEN -eq 1 ]; then
        echo "one genome, check quality"
        grep -v 'contamination' ${quality_csv} | tr ',' '\\t' | awk '{if(\$2 > 50 && \$3 < 5) print\$1}' > filtered_genomes.txt
        export COUNT=\$(less filtered_genomes.txt | wc -l)
        if [ \$COUNT -eq 1 ]; then
            cp genomes_folder/*.fa drep_output/dereplicated_genomes/
            ls drep_output/dereplicated_genomes > dereplicated_genomes.txt
        else
            echo "quality not passed -> no genomes"
        fi
    else
        dRep dereplicate -g genomes_folder/*.fa \
        -p ${task.cpus} \
        ${drep_params} \
        --genomeInfo ${quality_csv} \
        drep_output

        ls drep_output/dereplicated_genomes > dereplicated_genomes.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dRep: \$(echo \$(dRep -h 2>&1) | grep -o "dRep v[0-9.]\\+" | sed "s/dRep v//" )
    END_VERSIONS

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}
        genomes_folder: \$(ls genomes_folder | wc -l), dereplicated: \$(ls drep_output/dereplicated_genomes | wc -l)
    END_LOGGING
    """
}