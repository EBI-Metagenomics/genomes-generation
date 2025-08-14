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


/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_UNDERSCORE_TO_DOT {

    label 'process_low'
    tag "${contigs}"

    input:
    path(contigs)

    output:
    path("${contigs}"), emit: return_files

    script:
    """
    sed -i 's/\\_/\\./' ${contigs}
    """
}


process CHECKM2_TABLE_FOR_DREP_GENOMES {

    input:
    path(checkm_filtered_genomes_dir)
    path(dereplicated_genomes_tsv)

    output:
    path("checkm_results_mags.tab"), emit: checkm_results_mags

    script:
    """
    grep -f ${dereplicated_genomes_tsv} ${checkm_filtered_genomes_dir} > checkm_results_mags.tab || true
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
