/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT           } from '../../modules/local/align_bwa/main'
include { SAMTOOLS_FASTQ      } from '../../modules/nf-core/samtools/fastq/main'

workflow DECONTAMINATION {
    take:
        reads           // tuple(meta, contigs)
        ref_genome      // tuple(meta, ref_genome, ref_genome_index)

    main:

    ALIGNMENT(reads, ref_genome, false)

    SAMTOOLS_FASTQ(ALIGNMENT.out.bams.map{it -> [it[0], it[2]]}, false)

    emit:
        decontaminated_reads = SAMTOOLS_FASTQ.out.fastq
}