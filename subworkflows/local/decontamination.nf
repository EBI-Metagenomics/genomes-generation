/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT       } from '../../modules/local/align_bwa/main'
include { SAMTOOLS_BAM2FQ } from '../../modules/nf-core/samtools/bam2fq/main'

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

    SAMTOOLS_BAM2FQ( ALIGNMENT.out.bam.map { meta, ref_fasta, bam, bai -> [ meta, bam ] }, reads.map { meta, reads -> meta.single_end == false } )

    ch_versions = ch_versions.mix(ALIGNMENT.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    emit:
    decontaminated_reads = SAMTOOLS_BAM2FQ.out.reads
    versions = ch_versions                          // channel: [ versions.yml ]
}