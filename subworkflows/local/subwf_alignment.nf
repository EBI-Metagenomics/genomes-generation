/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { INDEX_FASTA_META      } from '../../modules/local/align_bwa'
include { ALIGNMENT             } from '../../modules/local/align_bwa'

workflow ALIGN {
    take:
        input_data  // tuple (meta, fasta, [reads])
    main:

    INDEX_FASTA_META(input_data.map(item -> tuple(item[0], item[1])))

    ALIGNMENT(
        input_data.map(item -> tuple(item[0], item[2])),
        INDEX_FASTA_META.out.fasta_with_index,
        true
    )

    emit:
        output = ALIGNMENT.out.bams  // [meta, fasta, bam, bai]
}