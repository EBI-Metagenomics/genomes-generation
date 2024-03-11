/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { INDEX_FASTA } from '../../modules/local/align_bwa/main'
include { ALIGNMENT_BAM   } from '../../modules/local/align_bwa/main'

workflow ALIGN {

    take:
    assembly_and_reads // tuple (meta, assembly_fasta, [reads])

    main:

    ch_versions = Channel.empty()

    INDEX_FASTA( assembly_and_reads.map { meta, assembly, reads -> [meta, assembly]} )

    reads_assembly_index = assembly_and_reads \
        .map { meta, assembly, reads -> [ meta, reads ] } \
        .join( INDEX_FASTA.out.fasta_with_index )

    ALIGNMENT_BAM(
        reads_assembly_index,
        true,
    )

    ch_versions = ch_versions.mix( INDEX_FASTA.out.versions.first() )
    ch_versions = ch_versions.mix( ALIGNMENT_BAM.out.versions.first() )

    emit:
    assembly_bam = ALIGNMENT_BAM.out.bam // [meta, assembly_fasta, bam, bai]
    versions     = ch_versions           // channel: [ versions.yml ]
}