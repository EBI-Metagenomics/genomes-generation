/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT_DECONTAMINATION    } from '../../modules/local/align_bwa/main'

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

    ALIGNMENT_DECONTAMINATION( to_align )

    ch_versions = ch_versions.mix(ALIGNMENT_DECONTAMINATION.out.versions.first())

    emit:
    decontaminated_reads = ALIGNMENT_DECONTAMINATION.out.reads
    versions = ch_versions                          // channel: [ versions.yml ]
}