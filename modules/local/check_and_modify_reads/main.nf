process CHECK_AND_MODIFY_READS {

    label 'process_medium'
    tag "${meta.id}"

    container 'quay.io/biocontainers/seqkit:2.7.0--h9ee0642_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*sanitized*.*.gz"), emit: modified_reads
    path "versions.yml"                      , emit: versions

    script:
    input_ch = reads.collect()
    if (meta.single_end ) {
        """
        seqkit sana ${input_ch[0]} | seqkit replace -p ${meta.id} -r ${meta.assembly_accession} -o ${meta.id}_sanitized.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version 2>&1 | sed 's/seqkit //g')
        END_VERSIONS
        """
    }
    else {
        """
        echo "read1"
        seqkit sana ${input_ch[0]} | seqkit replace -p ${meta.id} -r ${meta.assembly_accession} -o ${meta.id}_sanitized_1.fastq.gz

        echo "read2"
        seqkit sana ${input_ch[1]} | seqkit replace -p ${meta.id} -r ${meta.assembly_accession} -o ${meta.id}_sanitized_2.fastq.gz

        echo "sync read1 and read2"
        seqkit pair -1 ${meta.id}_sanitized_1.fastq.gz -2 ${meta.id}_sanitized_2.fastq.gz -O result

        mv result/${meta.id}_sanitized_1.fastq.gz ${meta.id}_sanitized_1.fastq.gz
        mv result/${meta.id}_sanitized_2.fastq.gz ${meta.id}_sanitized_2.fastq.gz
        rm -rf result

        echo "Done"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(seqkit version 2>&1 | sed 's/seqkit //g')
        END_VERSIONS
        """
    }
}