/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT      } from '../../modules/local/align_bwa/main'
include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'

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

    ALIGNMENT( to_align, false )

    SAMTOOLS_FASTQ( ALIGNMENT.out.bam.map{ it -> tuple( it[0], it[2] ) }, false )

    ch_versions = ch_versions.mix(ALIGNMENT.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first())

    emit:
    decontaminated_reads = SAMTOOLS_FASTQ.out.fastq // TODO: update
    versions = ch_versions                          // channel: [ versions.yml ]
}