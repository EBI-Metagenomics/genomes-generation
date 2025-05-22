/*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Decontamination subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT_DECONTAMINATION    } from '../../modules/local/bwa_mem2/decontamination/main'


workflow DECONTAMINATION {
    take:
    reads            // tuple(meta, reads)
    ref_genome       // path(reference_genome)
    ref_genome_index // path(reference_genome_index

    main:

    ch_versions = Channel.empty()
    if ( params.skip_decontamination ) {
        println('skipping decontamination')
        decontaminated_reads = reads
    }
    else {
        /*
        * decontaminate with bwa-mem2
        */
        ALIGNMENT_DECONTAMINATION(
           reads.map { meta, reads_item -> [ meta, reads_item, ref_genome, ref_genome_index ] }
        )
        ch_versions = ch_versions.mix(ALIGNMENT_DECONTAMINATION.out.versions.first())
        decontaminated_reads = ALIGNMENT_DECONTAMINATION.out.reads
    }

    emit:
    decontaminated_reads = decontaminated_reads
    versions             = ch_versions
}