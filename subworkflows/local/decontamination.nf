/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT_READS    } from '../../modules/local/align_bwa/main'

workflow DECONTAMINATION {
    take:
    reads            // tuple(meta, reads)
    ref_genome       // path(reference_genome)
    ref_genome_index // path(reference_genome_index

    main:

    ch_versions = Channel.empty()

    to_align = reads.map { meta, reads -> 
        [ meta, reads, ref_genome, ref_genome_index ]
    }

    ALIGNMENT_READS( to_align, false )

    ch_versions = ch_versions.mix(ALIGNMENT_READS.out.versions.first())

    emit:
    decontaminated_reads = ALIGNMENT_READS.out.reads
    versions = ch_versions                          // channel: [ versions.yml ]
}