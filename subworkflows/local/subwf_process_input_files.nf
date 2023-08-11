include { CHANGE_DOT_TO_UNDERSCORE_CONTIGS              } from '../../modules/local/utils'
include { CHANGE_DOT_TO_UNDERSCORE_READS                } from '../../modules/local/utils'
include { CHANGE_ERR_TO_ERZ_READS                       } from '../../modules/local/utils'


workflow PROCESS_INPUT {
    take:
        input_data            // tuple( meta, assembly_file, [raw_reads] )
        rename_file
    main:
    reads = input_data.map(item -> tuple(item[0], item[2]))
    contigs = input_data.map(item -> tuple(item[0], item[1]))

    // --- MODIFY CONTIGS
    // change . to _
    CHANGE_DOT_TO_UNDERSCORE_CONTIGS(contigs)           // tuple(meta, assembly)

    // --- MODIFY READS
    // change ERR in reads to ERZ
    CHANGE_ERR_TO_ERZ_READS(reads, rename_file.first())    // tuple(meta, [reads]])
    // change . to _
    CHANGE_DOT_TO_UNDERSCORE_READS(CHANGE_ERR_TO_ERZ_READS.out.return_files)

    emit:
        // tuple( meta, assembly_file, [raw_reads] )
        return_tuple = CHANGE_DOT_TO_UNDERSCORE_CONTIGS.out.underscore_contigs.combine(CHANGE_DOT_TO_UNDERSCORE_READS.out.underscore_reads, by: 0)
}
