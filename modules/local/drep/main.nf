/*
 * dRep, this workflow dereplicates a set of genomes.
*/
process DREP {

    tag "${meta.id} ${drep_params} ${type}"

    //conda "bioconda::drep=3.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/drep:3.2.2--pyhdfd78af_0':
        'quay.io/biocontainers/drep:3.2.2--pyhdfd78af_0' }"

    // that publishDir needs to be a function because of type parameter
    publishDir(
            path: "${params.outdir}/genomes_drep/${type}",
            mode: params.publish_dir_mode,
            failOnError: true,
            pattern: "dereplicated_genomes.txt",
            when: { ${meta.id} == "aggregated" },
    )
    publishDir(
            path: "${params.outdir}/genomes_drep",
            mode: params.publish_dir_mode,
            failOnError: true,
            pattern: "drep_output/dereplicated_genomes/*.fa",
            when: { ${meta.id} == "aggregated" },
            saveAs: { filename -> {
                        def output_file = new File(filename);
                        return "${type}/genomes/${output_file.name}";
                        }
            }
    )

    input:
    tuple val(meta), path(genomes_list, stageAs: "genomes_dir/*"), path(quality_csv)
    val drep_params
    val type

    output:
    tuple val(meta), path("drep_output/dereplicated_genomes/*"), optional: true, emit: dereplicated_genomes
    tuple val(meta), path("dereplicated_genomes.txt")          , optional: true, emit: dereplicated_genomes_list
    path "versions.yml"                                                        , emit: versions

    script:
    """
    mkdir -p drep_output/dereplicated_genomes
    LEN=\$(ls genomes_dir/*.fa | wc -l)
    if [ \$LEN -eq 0 ]; then
        echo "no genomes"
    elif [ \$LEN -eq 1 ]; then
        echo "one genome, check quality"
        grep -v 'contamination' ${quality_csv} | tr ',' '\\t' | awk '{if(\$2 > 50 && \$3 < 5) print\$1}' > filtered_genomes.txt
        COUNT=\$(less filtered_genomes.txt | wc -l)
        if [ \$COUNT -eq 1 ]; then
            cp genomes_dir/*.fa drep_output/dereplicated_genomes/
            ls drep_output/dereplicated_genomes > dereplicated_genomes.txt
        else
            echo "quality not passed -> no genomes"
        fi
    else
        dRep dereplicate -g genomes_dir/*.fa \
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
    """
}