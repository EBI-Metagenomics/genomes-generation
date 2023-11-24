process METABAT2_METABAT2 {
    tag "$meta.id ${fasta} with ${depth}"
    label 'process_medium'

    conda "bioconda::metabat2=2.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabat2:2.15--h986a166_1' :
        'biocontainers/metabat2:2.15--h986a166_1' }"

    input:
    tuple val(meta), path(fasta), path(depth)

    output:
    tuple val(meta), path("bins/*.tooShort.fa")                    , optional:true , emit: tooshort
    tuple val(meta), path("bins/*.lowDepth.fa")                    , optional:true , emit: lowdepth
    tuple val(meta), path("bins/*.unbinned.fa")                    , optional:true , emit: unbinned
    tuple val(meta), path("*.tsv.gz")                              , optional:true , emit: membership
    tuple val(meta), path("${meta.id}_metabat_bins")               , optional:true , emit: fasta
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def decompress_depth = depth ? "gzip -d -f $depth" : ""
    def depth_file = depth ? "-a ${depth.baseName}" : ""
    """
    $decompress_depth

    rm -rf metabat2
    rm -rf bins

    metabat2 \\
        $args \\
        -i $fasta \\
        $depth_file \\
        -t $task.cpus \\
        --saveCls \\
        -o metabat2/${prefix}

    mv metabat2/${prefix} ${prefix}.tsv
    mv metabat2 bins

    rm -f ${prefix}.tsv.gz
    gzip ${prefix}.tsv

    mkdir -p ${meta.id}_metabat_bins
    # I had to hardcode the version, it was failling with erro 141 for some reason
    version="2.15"

    if [ -z "\$(ls -A bins/*.[0-9]*.fa)" ]; then
        echo "Folder is empty"
    else
        for i in bins/*.[0-9]*.fa; do
            original_name=\$(basename \${i})
            new_filename=\${original_name#*-}
            digit=\${new_filename#*.}
            digit_without_extension=\${digit%.*}
            new_name="${meta.id}_metabat_bins/${meta.id}_bin.\${digit_without_extension}-metabat2_v\${version}.fa"
            mv \${i} \${new_name}
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS
    """
}