//
// Remove unwanted characters from the contigs
//

include { CHANGE_DOT_TO_UNDERSCORE_CONTIGS } from '../../modules/local/utils'
include { ERZ_TO_ERR                       } from '../../modules/local/utils'

workflow PROCESS_INPUT {

    take:
    input_data  // tuple( meta, assembly_file, [raw_reads] )

    main:

    ch_versions = Channel.empty()

    reads = input_data.map { meta, _, runs -> [meta, runs] }
    contigs = input_data.map { meta, assembly, _ -> [meta, assembly] }

    // --- MODIFY CONTIGS
    // change . to _
    CHANGE_DOT_TO_UNDERSCORE_CONTIGS( contigs ) // tuple(meta, assembly)

    // --- MODIFY READS
    // change ERR in reads to ERZ
    // TODO: make this process faster (biopython is slow)
    ERZ_TO_ERR( reads ) // input: tuple(meta, [reads]]); output: tuple(meta, [modified_reads]])

    result = CHANGE_DOT_TO_UNDERSCORE_CONTIGS.out.underscore_contigs.join(
        ERZ_TO_ERR.out.modified_reads
    )

    ch_versions.mix( ERZ_TO_ERR.out.versions.first() )

    emit:
    assembly_and_reads = result // tuple( meta, assembly_file, [raw_reads] )
    versions           = ch_versions
}
