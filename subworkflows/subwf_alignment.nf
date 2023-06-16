/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT_WITH_INDEXING } from '../modules/align_bwa'

workflow ALIGN {
    take:
        input_data  // tuple (name, contig, [reads])

    main:

    sample_reads = input_data.map(item -> tuple(item[0], item[2]))
    ref_db = input_data.map(item -> item[1])
    getFastaBasename = { fasta_file ->
            def bname = fasta_file[1].toString().tokenize("/")[-1]
            return bname
        }
    assembly = input_data.map(getFastaBasename)
    samtools_args = channel.value("-q 20 -Sb")
    ALIGNMENT_WITH_INDEXING(sample_reads, ref_db, assembly, samtools_args)

    emit:
        annotated_bams = ALIGNMENT_WITH_INDEXING.out.bams
}