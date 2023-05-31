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
    input:
    path(file_to_compress)

    output:
    path("*.gz"), emit: compressed

    script:
    """
    gzip ${file_to_compress}
    """
}

/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_DOT_TO_UNDERSCORE {

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    publishDir(
        path: "${params.outdir}/prepare_data",
        mode: 'copy',
        failOnError: true
    )

    input:
    tuple val(accession), path(contigs)

    output:
    tuple val(accession), path("${contigs}"), emit: return_contigs

    script:
    """
    sed -i 's/\\./\\_/' ${contigs}
    """
}

/*
 * clean reads, change . to _ from contigs
*/
process CHANGE_UNDERSCORE_TO_DOT {

    input:
    path contigs

    output:
    path "${contigs.baseName}", emit: return_contigs

    script:
    """
    sed -i 's/\\_/\\./' ${contigs}
    """
}

/*
 * change run accession to assembly accession
*/
process CHANGE_ERR_TO_ERZ {

    input:
    tuple val(run_accession), path(files)
    path rename_file

    output:
    tuple val(run_accession), path("*changed*.fastq.gz"), emit: return_files

    script:
    input_ch = files.collect()
    if (input_ch.size() == 1 ) {
        print('single')
        """
        echo ${run_accession}
        grep "${run_accession}" ${rename_file} > help_file
        export from_accession=\$(cat help_file | cut -f1)
        export to_accession=\$(cat help_file | cut -f2)
        echo "\${from_accession} --> \${to_accession}"

        zcat "${input_ch[0]}" | sed "s/\${from_accession}/\${to_accession}/g" | gzip > ${run_accession}_changed.fastq.gz
        """
    }
    else if (input_ch.size() == 2 ) {
        """
        echo ${run_accession}
        grep "${run_accession}" ${rename_file} > help_file
        export from_accession=\$(cat help_file | cut -f1)
        export to_accession=\$(cat help_file | cut -f2)
        echo "\${from_accession} --> \${to_accession}"

        zcat "${input_ch[0]}" | sed "s/\${from_accession}/\${to_accession}/g" | gzip > ${run_accession}_changed_1.fastq.gz
        zcat "${input_ch[1]}" | sed "s/\${from_accession}/\${to_accession}/g" | gzip > ${run_accession}_changed_2.fastq.gz
        """
    }
    else {
        print('incorrect input') }
}