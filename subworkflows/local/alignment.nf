/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { INDEX_FASTA }     from '../../modules/local/align_bwa/main'
include { FEATURED_ALIGNMENT   } from '../../modules/local/align_bwa/main'

workflow ALIGN {

    take:
    assembly_and_reads // tuple (meta, assembly_fasta, [reads])
    run_jgi_depth
    run_concoct_table_generation
    concoct_generate_bed_file

    main:

    ch_versions = Channel.empty()

    INDEX_FASTA( assembly_and_reads.map { meta, assembly, reads -> [meta, assembly]} )

    reads_assembly_index = assembly_and_reads \
        .map { meta, assembly, reads -> [ meta, reads ] } \
        .join( INDEX_FASTA.out.fasta_with_index )

    FEATURED_ALIGNMENT( reads_assembly_index, run_jgi_depth, run_concoct_table_generation, concoct_generate_bed_file )

    ch_versions = ch_versions.mix( INDEX_FASTA.out.versions.first() )
    ch_versions = ch_versions.mix( FEATURED_ALIGNMENT.out.versions.first() )

    emit:
    jgi_depth         = FEATURED_ALIGNMENT.out.depth  // [ meta, depth.txt.gz ]
    concoct_data      = FEATURED_ALIGNMENT.out.concoct_data  // [ meta, concoct.tsv ]
    samtools_idxstats = FEATURED_ALIGNMENT.out.idxstats
    versions          = ch_versions           // channel: [ versions.yml ]
}