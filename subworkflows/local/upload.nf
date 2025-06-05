//
// Create table for genome uploader and upload mags
//

include { CREATE_MANIFESTS_FOR_UPLOAD } from '../../modules/local/genome_uploader/create_manifests'
include { CHECK_WEBIN_SUCCESS         } from '../../modules/local/genome_uploader/check_webin_success'
include { PREPARE_TSV_FOR_UPLOADER    } from '../../modules/local/genome_uploader/generate_table_for_upload/main'
include { WEBIN_CLI_UPLOAD            } from '../../modules/local/genome_uploader/webin_cli_upload'


workflow UPLOAD_MAGS {

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

    if (params.upload_mags || params.upload_bins) {
        mags = euk_genomes.combine(prok_genomes).filter { it != null }
        CREATE_MANIFESTS_FOR_UPLOAD(
           PREPARE_TSV_FOR_UPLOADER.out.tsv_for_uploader,
           mags
        )
        // combine mag with corresponding manifest
        mags_ch = mags.flatten().map { bin ->
           def prefix = bin.name.replaceAll(/\.fa\.gz$/, '')
           [prefix, bin]
        }
        manifests_ch = CREATE_MANIFESTS_FOR_UPLOAD.out.manifests.flatten()
         .map { manifest ->
            def prefix = manifest.name.replaceAll(/_\d+\.manifest$/, '')
            [prefix, manifest]
         }
        combined_ch = mags_ch.join(manifests_ch)

        WEBIN_CLI_UPLOAD(combined_ch)
        ch_versions = ch_versions.mix( WEBIN_CLI_UPLOAD.out.versions )

        CHECK_WEBIN_SUCCESS(WEBIN_CLI_UPLOAD.out.webin_report)
        // Count successes using channel operators
        success_count = CHECK_WEBIN_SUCCESS.out.results
           .map { sample_id, status -> status == "true" ? 1 : 0 }
           .sum()
           .view { count -> "✅ Total successful submissions: ${count}" }

        // Count total submissions
        total_count = CHECK_WEBIN_SUCCESS.out.results
           .count()
           .view { count -> "📊 Total submissions processed: ${count}" }

    }

    emit:
    versions = ch_versions
}