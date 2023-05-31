/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { INDEX_REF_GENOME } from '../modules/align_bwa'
include { ALIGNMENT } from '../modules/align_bwa'

workflow ALIGN {
    take:
        contigs  // tuple(name, contigs)
        reads    // tuple(name, contigs)

    main:
    INDEX_REF_GENOME(contigs)

    ref_folders = INDEX_REF_GENOME.out.index_folder
    sample_data = reads.combine(ref_folders, by: 0).combine(contigs, by: 0)  // because it can swap ref_db and reads
    sample_reads = sample_data.map(item -> tuple(item[0], item[1]))
    ref_db = sample_data.map(item -> item[2])
    getFastaBasename = { fasta_file ->
            def bname = fasta_file[3].toString().tokenize("/")[-1]
            return bname
        }
    assembly = sample_data.map(getFastaBasename)
    ALIGNMENT(sample_reads, ref_db, assembly)

    emit:
        annotated_bams = ALIGNMENT.out.bams
}