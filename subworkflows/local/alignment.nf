/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { INDEX_FASTA } from '../../modules/local/align_bwa/main'
include { ALIGNMENT   } from '../../modules/local/align_bwa/main'

workflow ALIGN {

    take:
    assembly_and_reads  // tuple (meta, assembly_fasta, [reads])

    main:

    INDEX_FASTA( assembly_and_reads.map(item -> tuple(item[0], item[1])) )

    to_align = assembly_and_reads \
        .map(item -> tuple(item[0], item[2])) \
        .join(INDEX_FASTA.out.fasta_with_index, by: [0]) 

    ALIGNMENT(
        to_align,
        true,
    )

    ch_versions = ch_versions.mix(INDEX_FASTA.out.versions.first())
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions.first())

    emit:
    output = ALIGNMENT.out.bams  // [meta, assembly_fasta, bam, bai]
    versions = ch_versions       // channel: [ versions.yml ]
}