//
// Create table for genome uploader
//

include { PREPARE_TSV_FOR_UPLOADER    } from '../../modules/local/genome_uploader/generate_table_for_upload/main'

workflow PREPARE_UPLOAD_FILES {

    take:
    euk_genomes
    prok_genomes
    assembly_software
    stats_euks
    stats_proks
    coverage_euks
    coverage_proks
    rna
    taxonomy_euks
    taxonomy_proks

    main:

    ch_versions = Channel.empty()

    PREPARE_TSV_FOR_UPLOADER(
            euk_genomes,
            prok_genomes,
            assembly_software,
            stats_euks,
            stats_proks,
            coverage_euks,
            coverage_proks,
            rna,
            taxonomy_euks,
            taxonomy_proks
        )
    ch_versions = ch_versions.mix( PREPARE_TSV_FOR_UPLOADER.out.versions )

    emit:
    versions = ch_versions
}