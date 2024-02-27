//
// Remove unwanted characters from the contigs
//

include { CHANGE_DOT_TO_UNDERSCORE_CONTIGS } from '../../modules/local/utils'
include { ERR_TO_ERZ                       } from '../../modules/local/utils'

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
    ERR_TO_ERZ( reads ) // input: tuple(meta, [reads]]); output: tuple(meta, [modified_pe_reads]/modified_se_reads])

    // this process returns list of reads if they are PE and path for SE
    // it is required for next step FASTP that SE reads should be also in list
    reads_changed = ERR_TO_ERZ.out.modified_reads.map{meta, reads ->
                                                            result = [meta]
                                                            if (reads.collect().size() == 2) {
                                                                result.add(reads) }
                                                            else {
                                                                result.add([reads]) }
                                                            return result
                                                          }

    result = CHANGE_DOT_TO_UNDERSCORE_CONTIGS.out.underscore_contigs.join(reads_changed)

    ch_versions = ch_versions.mix( ERR_TO_ERZ.out.versions.first() )

    emit:
    assembly_and_reads = result // tuple( meta, assembly_file, [raw_reads] )
    versions           = ch_versions
}
