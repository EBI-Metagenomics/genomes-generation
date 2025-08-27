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
    mags_or_bins_flag

    main:

    ch_versions = Channel.empty()
    
    /* --  Generate TSV for genome_uploader -- */
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

    mags = euk_genomes.combine(prok_genomes).filter { it != null }

    /* --  Generate manifests for submission -- */
    CREATE_MANIFESTS_FOR_UPLOAD(
        PREPARE_TSV_FOR_UPLOADER.out.tsv_for_uploader,
        mags,
        mags_or_bins_flag
    )
    ch_versions = ch_versions.mix( CREATE_MANIFESTS_FOR_UPLOAD.out.versions )

    /* --   combine MAG with corresponding manifest -- */
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

    /* --   Upload MAGs with webin-cli -- */
    WEBIN_CLI_UPLOAD(combined_ch)
    ch_versions = ch_versions.mix( WEBIN_CLI_UPLOAD.out.versions )

    // Count successes using channel operators
    WEBIN_CLI_UPLOAD.out.upload_status
        .map { _sample_id, status -> status == "success" ? 1 : 0 }
        .sum()
        .view { count -> "✅ Total successful submissions: ${count}" }

    // Count total submissions
    WEBIN_CLI_UPLOAD.out.upload_status
        .count()
        .view { count -> "📊 Total submissions processed: ${count}" }

    // Generate summary report file
    WEBIN_CLI_UPLOAD.out.upload_status
        .map { _sample_id, status -> status == "success" ? 1 : 0 }
        .sum()
        .combine(
            WEBIN_CLI_UPLOAD.out.upload_status.count()
        )
        .map { success_count_value, total_count_value ->
            """ENA Submission Summary Report
    Generated: ${new Date()}
    Total successful submissions: ${success_count_value}
    Total submissions processed: ${total_count_value}
    Status: ${success_count_value == total_count_value ? 'ALL SUBMISSIONS SUCCESSFUL' : 'SOME SUBMISSIONS FAILED'}
    """
        }
        .collectFile(name: 'ena_submission_summary.txt', storeDir: "${params.outdir}/upload")

    emit:
    versions = ch_versions
}