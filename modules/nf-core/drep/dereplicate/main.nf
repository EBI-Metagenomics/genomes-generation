process DREP_DEREPLICATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/drep:3.6.2--pyhdfd78af_0'
        : 'biocontainers/drep:3.6.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fastas, stageAs: 'input_fastas/*'), path(quality_csv)
    tuple val(meta2), path(drep_work, stageAs: 'drep_work/')

    output:
    tuple val(meta), path("dereplicated_genomes/*"), optional: true, emit: fastas
    tuple val(meta), path("data_tables/*.csv")     , optional: true, emit: summary_tables
    tuple val(meta), path("figures/*pdf")          , optional: true, emit: figures
    tuple val(meta), path("logger.log")            , optional: true, emit: log
    path "versions.yml"                                            , emit: versions
    path "progress.log"                                            , emit: progress_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ ! -d drep_work/ ]]; then
        mkdir drep_work/
    fi

    find -L input_fastas/ -type f > fastas_paths.txt
    GENOME_COUNT=\$(cat fastas_paths.txt | wc -l)

    if [[ \$GENOME_COUNT -eq 0 ]]; then
        echo "No genomes found in input_fastas/"

    elif [[ \$GENOME_COUNT -eq 1 ]]; then
        echo "Only one genome detected. Skipping dRep and copying genome directly."
        
        # Copy the single genome to dereplicated_genomes
        cp input_fastas/* dereplicated_genomes/
        
        # Create minimal output files required by PROPAGATE_TAXONOMY_TO_BINS process
        GENOME_NAME=\$(basename input_fastas/*)
        echo "genome,cluster" > data_tables/Wdb.csv
        echo "\$GENOME_NAME,1_0" >> data_tables/Wdb.csv
        echo "genome,secondary_cluster" > data_tables/Cdb.csv
        echo "\$GENOME_NAME,1_0" >> data_tables/Cdb.csv
        
    else
        echo "Multiple genomes detected (\$GENOME_COUNT). Running dRep."
        
        dRep \\
            dereplicate \\
            drep_work/ \\
            -p ${task.cpus} \\
            -g fastas_paths.txt \\
            --genomeInfo ${quality_csv} \
            ${args} \\

        ## We copy the output files to copies to ensure we don't break an already
        ## existing drep_work/ from an upstream dRep module (e.g. compare)
        mkdir dereplicated_genomes/ figures/ data_tables
        cp drep_work/dereplicated_genomes/* dereplicated_genomes/
        cp drep_work/figures/* figures/
        cp drep_work/data_tables/* data_tables/
        cp drep_work/log/logger.log logger.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(dRep | head -n 2 | sed 's/.*v//g;s/ .*//g' | tail -n 1)
    END_VERSIONS

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}
        genomes_folder: \$(ls genomes_folder | wc -l), dereplicated: \$(ls drep_output/dereplicated_genomes | wc -l)
    END_LOGGING
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir ${prefix}/
    mkdir -p dereplicated_genomes/ figures/ data_tables/
    touch dereplicated_genomes/{test1,test2}.fasta
    touch figures/{Clustering_scatterplots,Cluster_scoring,Primary_clustering_dendrogram,Secondary_clustering_dendrograms,Winning_genomes}.pdf
    touch data_tables/{Bdb,Cdb,Chdb,Mdb,Ndb,Sdb,Wdb,Widb}.csv
    touch logger.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(dRep | head -n 2 | sed 's/.*v//g;s/ .*//g' | tail -n 1)
    END_VERSIONS

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}
        input_fastas: \$(ls input_fastas | wc -l), dereplicated: \$(ls dereplicated_genomes/ | wc -l)
    END_LOGGING
    """
}
