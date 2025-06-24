//
// Create table for genome uploader and upload mags
//

include { CREATE_MANIFESTS_FOR_UPLOAD } from '../../modules/local/genome_uploader/create_manifests'
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
            def prefix = manifest.name.replaceAll(/\.manifest$/, '')
            [prefix, manifest]
         }
        combined_ch = mags_ch.join(manifests_ch)

        WEBIN_CLI_UPLOAD(combined_ch)
        ch_versions = ch_versions.mix( WEBIN_CLI_UPLOAD.out.versions )

        // Count successes using channel operators
        success_count = WEBIN_CLI_UPLOAD.out.upload_status
           .map { sample_id, status -> status == "true" ? 1 : 0 }
           .sum()
           .view { count -> "✅ Total successful submissions: ${count}" }

        // Count total submissions
        total_count = WEBIN_CLI_UPLOAD.out.upload_status
           .count()
           .view { count -> "📊 Total submissions processed: ${count}" }

        // Generate summary report file
        WEBIN_CLI_UPLOAD.out.upload_status
            .map { sample_id, status -> status == "true" ? 1 : 0 }
            .sum()
            .combine(
                WEBIN_CLI_UPLOAD.out.upload_status.count()
            )
            .map { success_count, total_count ->
                """ENA Submission Summary Report
        Generated: ${new Date()}
        Total successful submissions: ${success_count}
        Total submissions processed: ${total_count}
        Status: ${success_count == total_count ? 'ALL SUBMISSIONS SUCCESSFUL' : 'SOME SUBMISSIONS FAILED'}
        """
            }
            .collectFile(name: 'ena_submission_summary.txt', storeDir: "${params.outdir}/upload")
            }

    emit:
    versions = ch_versions
}
