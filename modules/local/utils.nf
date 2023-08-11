/*
    ~~~~~~~~~~~~~~~~~~
     uncompress
    ~~~~~~~~~~~~~~~~~~
*/
process GUNZIP {
    input:
    path(compressed_file)

    output:
    path("out/*"), emit: uncompressed

    script:
    """
    mkdir out
    cp ${compressed_file} out
    cd out
    gunzip *
    """
}

/*
    ~~~~~~~~~~~~~~~~~~
     compress
    ~~~~~~~~~~~~~~~~~~
*/
process GZIP {
    stageInMode 'copy'
    tag "${file_to_compress}"

    publishDir(
        path: "${params.outdir}/${output_folder}",
        mode: 'copy',
        failOnError: true
    )

    input:
    file(file_to_compress)
    val output_folder

    output:
    path("*.gz"), emit: compressed

    script:
    """
    pigz ${file_to_compress}
    """
}

/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_DOT_TO_UNDERSCORE_CONTIGS {

    tag "${meta.id}"
    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    publishDir(
        path: "${params.outdir}/intermediate_steps/prepare_data",
        mode: 'copy',
        failOnError: true
    )

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${contigs}"), emit: underscore_contigs

    script:
    """
    sed -i 's/\\./\\_/' ${contigs}
    """
}

process CHANGE_DOT_TO_UNDERSCORE_READS {
    tag "${meta.id}"
    container 'quay.io/biocontainers/pigz:2.3.4'

    publishDir(
        path: "${params.outdir}/intermediate_steps/prepare_data",
        mode: 'copy',
        failOnError: true
    )

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*underscore*.fastq.gz"), emit: underscore_reads

    script:
    input_ch = reads.collect()
    if (input_ch.size() == 1 ) {
        """
        zcat "${input_ch[0]}" | awk '{if (NR%4==1){gsub(/\\./,"_")}; print}' | pigz > ${meta.id}_underscore.fastq.gz
        """
    }
    else if (input_ch.size() == 2 ) {
        """
        zcat "${input_ch[0]}" | awk '{if (NR%4==1){gsub(/\\./,"_")}; print}' | pigz > ${meta.id}_underscore_1.fastq.gz
        zcat "${input_ch[1]}" | awk '{if (NR%4==1){gsub(/\\./,"_")}; print}' | pigz > ${meta.id}_underscore_2.fastq.gz
        """
    }
    else {
        print('incorrect input') }
}

/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_UNDERSCORE_TO_DOT {

    tag "${name} ${contigs}"

    input:
    tuple val(name), path(contigs)

    output:
    tuple val(name), path("${contigs}"), emit: return_files

    script:
    """
    sed -i 's/\\_/\\./' ${contigs}
    """
}

/*
 * change run accession to assembly accession
*/
process CHANGE_ERR_TO_ERZ_READS {
    tag "${meta.id}"
    container 'quay.io/biocontainers/pigz:2.3.4'

    input:
    tuple val(meta), path(files)
    path rename_file

    output:
    tuple val(meta), path("*changed*.fastq.gz"), emit: return_files

    script:
    input_ch = files.collect()
    if (input_ch.size() == 1 ) {
        """
        echo ${meta.id}
        grep "${meta.id}" ${rename_file} > help_file
        export from_accession=\$(cat help_file | cut -f1)
        export to_accession=\$(cat help_file | cut -f2)
        echo "\${from_accession} --> \${to_accession}"

        zcat "${input_ch[0]}" | sed 's/\${from_accession}/\${to_accession}/g' | pigz > ${meta.id}_changed.fastq.gz
        """
    }
    else if (input_ch.size() == 2 ) {
        """
        echo ${meta.id}
        grep "${meta.id}" ${rename_file} > help_file
        export from_accession=\$(cat help_file | cut -f1)
        export to_accession=\$(cat help_file | cut -f2)
        echo "\${from_accession} --> \${to_accession}"

        zcat "${input_ch[0]}" | sed 's/\${from_accession}/\${to_accession}/g' | pigz > ${meta.id}_changed_1.fastq.gz
        zcat "${input_ch[1]}" | sed 's/\${from_accession}/\${to_accession}/g' | pigz > ${meta.id}_changed_2.fastq.gz
        """
    }
    else {
        print('incorrect input') }
}