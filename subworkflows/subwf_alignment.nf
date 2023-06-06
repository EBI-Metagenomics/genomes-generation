/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT_WITH_INDEXING } from '../modules/align_bwa'

workflow ALIGN {
    take:
        contigs  // tuple(name, contigs)
        reads    // tuple(name, contigs)

    main:

    sample_data = reads.combine(contigs, by: 0)  // because it can swap contigs and reads
    sample_reads = sample_data.map(item -> tuple(item[0], item[1]))
    ref_db = sample_data.map(item -> item[2])
    getFastaBasename = { fasta_file ->
            def bname = fasta_file[2].toString().tokenize("/")[-1]
            return bname
        }
    assembly = sample_data.map(getFastaBasename)
    samtools_args = channel.value("-q 20 -Sb")
    ALIGNMENT_WITH_INDEXING(sample_reads, ref_db, assembly, samtools_args)

    emit:
        annotated_bams = ALIGNMENT_WITH_INDEXING.out.bams
}