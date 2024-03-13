/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { INDEX_FASTA } from '../../modules/local/align_bwa/main'
include { ALIGNMENT_DEPTH   } from '../../modules/local/align_bwa/main'

workflow ALIGN {

    take:
    assembly_and_reads // tuple (meta, assembly_fasta, [reads])

    main:

    ch_versions = Channel.empty()

    INDEX_FASTA( assembly_and_reads.map { meta, assembly, reads -> [meta, assembly]} )

    reads = assembly_and_reads \
        .map { meta, assembly, reads -> [ meta, reads ] }
    reads_assembly_index = reads.join( INDEX_FASTA.out.fasta_with_index )

    ALIGNMENT_DEPTH( reads_assembly_index )

    ch_versions = ch_versions.mix( INDEX_FASTA.out.versions.first() )
    ch_versions = ch_versions.mix( ALIGNMENT_DEPTH.out.versions.first() )

    emit:
    assembly_depth_reads = ALIGNMENT_DEPTH.out.depth.join(reads) // [meta, assembly_fasta, depth, [reads]]
    versions     = ch_versions           // channel: [ versions.yml ]
}