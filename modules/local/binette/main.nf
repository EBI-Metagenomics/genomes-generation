process BINETTE {

    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_binette:9744102b28ea4ec5' :
        'community.wave.seqera.io/library/binette_checkm2_diamond_pyfastx_pruned:4e303ff036c20212' }"
    
    input:
    tuple val(meta), path(bin_dirs), path(assembly)

    output:
    path "", emit: all_bin2classification_human_readable
    path "" , emit: all_bin2classification
    path ""               , emit: versions

    script:
    """
    binette --bin_dirs $bin_dirs --contigs $assembly

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        binette: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
