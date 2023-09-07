//
// Process the reads with FastP
// - QC reads
// - merge pairs
//

include { FASTP as FASTP_SINGLE       } from '../../../modules/nf-core/fastp/main'
include { FASTP as FASTP_PAIRED       } from '../../../modules/nf-core/fastp/main'

workflow QC_AND_MERGE {

    // Taken from: https://github.com/nf-core/taxprofiler/blob/master/subworkflows/local/shortread_fastp.nf

    take:
    reads // [[meta], [reads]]

    main:
    ch_versions           = Channel.empty()
    ch_multiqc_files      = Channel.empty()

    ch_input_for_fastp = reads.branch{
        single: it[0].single_end == true
        paired: it[0].single_end == false
    }

    // We don't provide the adapter sequences, which is the second parameter for fastp

    FASTP_SINGLE ( ch_input_for_fastp.single, [], false, false )

    // Last parameter here turns on merging of PE data
    FASTP_PAIRED ( ch_input_for_fastp.paired, [], false, params.shortread_qc_mergepairs )

    if ( params.shortread_qc_mergepairs ) {

        ch_fastp_reads_prepped_pe = FASTP_PAIRED.out.reads_merged.map {
            meta, reads ->
                def meta_new = meta.clone()
                meta_new['single_end'] = true
                [ meta_new, [ reads ].flatten() ]
        }
        ch_fastp_reads_prepped = ch_fastp_reads_prepped_pe.mix(
            FASTP_SINGLE.out.reads
        )

    } else {
        ch_fastp_reads_prepped = FASTP_PAIRED.out.reads.mix(
            FASTP_SINGLE.out.reads
        )
    }

    ch_versions = ch_versions.mix(FASTP_SINGLE.out.versions.first())
    ch_versions = ch_versions.mix(FASTP_PAIRED.out.versions.first())

    ch_processed_reads = ch_fastp_reads_prepped

    ch_multiqc_files = ch_multiqc_files.mix( FASTP_SINGLE.out.json )
    ch_multiqc_files = ch_multiqc_files.mix( FASTP_PAIRED.out.json )

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files

}