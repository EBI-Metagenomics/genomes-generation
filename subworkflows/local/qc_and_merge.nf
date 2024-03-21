//
// Process the reads with FastP
// - QC reads
// - merge pairs
//

include { FASTP as FASTP_SINGLE } from '../../modules/nf-core/fastp/main'
include { FASTP as FASTP_PAIRED } from '../../modules/nf-core/fastp/main'

workflow QC_AND_MERGE_READS {

    // Taken from: https://github.com/nf-core/taxprofiler/blob/master/subworkflows/local/shortread_fastp.nf

    take:
    reads // [[meta], [reads]]

    main:
    ch_versions           = Channel.empty()
    ch_multiqc_files      = Channel.empty()

    reads.branch {
        single: it[0].single_end
        paired: !it[0].single_end
    }.set {
        ch_input_for_fastp
    }

    // We don't provide the adapter sequences, which is the second parameter for fastp

    FASTP_SINGLE ( ch_input_for_fastp.single, [], false, false )

    // Last parameter here turns on merging of PE data
    FASTP_PAIRED ( ch_input_for_fastp.paired, [], false, params.merge_pairs )

    if ( params.merge_pairs ) {
        ch_fastp_paired = FASTP_PAIRED.out.reads_merged
    } else {
        ch_fastp_paired = FASTP_PAIRED.out.reads
    }

    ch_processed_reads = ch_fastp_paired.mix(
            FASTP_SINGLE.out.reads
        )

    ch_versions = ch_versions.mix(FASTP_SINGLE.out.versions.first())
    ch_versions = ch_versions.mix(FASTP_PAIRED.out.versions.first())

    ch_multiqc_files = ch_multiqc_files.mix( FASTP_SINGLE.out.json )
    ch_multiqc_files = ch_multiqc_files.mix( FASTP_PAIRED.out.json )

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files

}