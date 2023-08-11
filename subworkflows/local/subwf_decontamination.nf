/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT           } from '../../modules/local/align_bwa'
include { SAMTOOLS_FASTQ      } from '../../modules/nf-core/samtools/fastq/main'
include { GZIP                } from '../../modules/local/utils'

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