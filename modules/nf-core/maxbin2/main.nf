process MAXBIN2 {
    tag "$meta.id $contigs"
    label 'process_medium'

    conda "bioconda::maxbin2=2.2.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/maxbin2:2.2.7--hdbdd923_5' :
        'biocontainers/maxbin2:2.2.7--he1b5a44_2' }"

    input:
    tuple val(meta), path(contigs), path(reads), path(abund)

    output:
    tuple val(meta), path("maxbin_output"), emit: binned_fastas
    tuple val(meta), path("*.summary")    , emit: summary,  optional: true
    tuple val(meta), path("*.log.gz")     , emit: log
    tuple val(meta), path("*.marker.gz")  , emit: marker_counts, optional: true
    tuple val(meta), path("*.noclass.gz") , emit: unbinned_fasta, optional: true
    tuple val(meta), path("*.tooshort.gz"), emit: tooshort_fasta, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def associate_files = reads ? "-reads $reads" : "-abund $abund"
    """

    mkdir -p maxbin_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        maxbin2: \$( run_MaxBin.pl -v | head -n 1 | sed 's/MaxBin //' )
    END_VERSIONS

    set +e
    run_MaxBin.pl \\
        -contig $contigs \\
        $associate_files \\
        -thread $task.cpus \\
        $args \\
        -out $prefix

    MAXBIN_EXITCODE="\$?"
    set -e

    echo "Exit code \$MAXBIN_EXITCODE"
    if [ "\$MAXBIN_EXITCODE" == "255" ]; then
        echo "Maxbin2 exit code \$MAXBIN_EXITCODE"
        if [ ! -e "${prefix}.log" ]; then
            echo "maxbin.log does not exist."
            if [ ! -e "${prefix}.tooshort" ]; then
                echo "Unknown error" >&2
                exit \$MAXBIN_EXITCODE
            else
                echo "Checking for all contigs were moved to .tooshort file"
                export INITIAL_NUMBER=\$(grep '>' $contigs | wc -l)
                export TOO_SHORT_NUMBER=\$(grep '>' ${prefix}.tooshort | wc -l)
                if [[ \$INITIAL_NUMBER -eq \$TOO_SHORT_NUMBER ]]; then
                    echo "All contigs are short"
                    gzip *log
                    exit 0
                else
                    echo "Not all contigs are short. Exit with unknown error"
                    exit \$MAXBIN_EXITCODE
                fi
            fi
        else
            echo "maxbin.log exists -> checking for not enough marker genes error"
            if grep -q "Marker gene search reveals that the dataset cannot be binned (the medium of marker gene number <= 1). Program stop." "${prefix}.log"; then
                echo "Not enough marker genes"
                gzip *log
                exit 0
            else
                echo "Unknown problem. Exit"
                exit \$MAXBIN_EXITCODE
            fi
        fi
    fi

    echo "Collect folder"
    mv $prefix*.fasta maxbin_output/  || true
    if [ -z "\$(ls -A maxbin_output)" ]; then
        echo "Folder is empty"
    else
        for i in maxbin_output/*; do
            number=\$(echo "\$i" | cut -d '.' -f 2 | sed 's/^0*//')
            new_name="maxbin_output/${meta.id}_maxbin2_\${number}.fa"
            echo "\${i} to \${new_name}"
            mv \${i} \${new_name}
        done
    fi
    echo "Compress files"
    gzip *.noclass *.tooshort *log *.marker || true
    """
}
