//
// Remove unwanted characters from the contigs
//

include { CHANGE_DOT_TO_UNDERSCORE_CONTIGS  } from '../../modules/local/utils'
include { ERZ_TO_ERR                        } from '../../modules/local/utils'

workflow PROCESS_INPUT {

    take:
    input_data  // tuple( meta, assembly_file, [raw_reads] )
    erz_to_err_mapping_file

    main:

    reads = input_data.map(item -> tuple(item[0], item[2]))
    contigs = input_data.map(item -> tuple(item[0], item[1]))

    // --- MODIFY CONTIGS
    // change . to _
    CHANGE_DOT_TO_UNDERSCORE_CONTIGS(contigs) // tuple(meta, assembly)

    // --- MODIFY READS
    // change ERR in reads to ERZ
    // TODO: make this process faster (biopython is slow)
    ERZ_TO_ERR( reads, erz_to_err_mapping_file ) // tuple(meta, [reads]])

    result = CHANGE_DOT_TO_UNDERSCORE_CONTIGS.out.underscore_contigs.combine(
        ERZ_TO_ERR.out.modified_reads, by: 0
    )

    emit:

    // tuple( meta, assembly_file, [raw_reads] )
    assembly_and_reads = result
}
