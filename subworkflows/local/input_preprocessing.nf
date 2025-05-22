//
// Change reads and assemblies input files if necessary
//

include { FASTQC         } from '../../modules/nf-core/fastqc/main'
include { MODIFY_CONTIGS } from '../../modules/local/utils'
include { ERR_TO_ERZ     } from '../../modules/local/utils'


workflow INPUT_PREPROCESSING {

    take:
    input_data  // tuple( meta, assembly_file, [raw_reads] )

    main:

    ch_versions = Channel.empty()

    reads = input_data.map { meta, _, runs -> [meta, runs] }
    contigs = input_data.map { meta, assembly, _ -> [meta, assembly] }

    //
    // --- check input reads quality ---
    //
    FASTQC(
        reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // --- skip reads and contigs modifications --
    //
    if ( params.skip_preprocessing_input ) {
        println('skipping pre-processing')
        result = input_data
    }
    //
    // --- add correct accessions
    //
    else {
        /*
        * --- MODIFY CONTIGS ( change . to _ )
        *  output: tuple(meta, assembly)
        */
        MODIFY_CONTIGS(
            contigs
        )

        /*
        * --- MODIFY READS ( change ERR in reads to ERZ )
        *  input: tuple(meta, [reads]])
        *  output: tuple(meta, [modified_pe_reads]/modified_se_reads])
        */
        ERR_TO_ERZ(
            reads
        )
        ch_versions = ch_versions.mix( ERR_TO_ERZ.out.versions.first() )

        // this process returns list of reads if they are PE and path for SE
        // it is required for next step FASTP that SE reads should be also in list
        reads_changed = ERR_TO_ERZ.out.modified_reads.map{ meta, reads ->
                                                                result = [meta]
                                                                if (!meta.single_end) {
                                                                    result.add(reads) }
                                                                else {
                                                                    result.add([reads]) }
                                                                return result
                                                              }

        result = MODIFY_CONTIGS.out.underscore_contigs.join(reads_changed)
    }

    emit:
    assembly_and_reads = result // tuple( meta, assembly_file, [raw_reads] )
    fastqc_gzip        = FASTQC.out.zip
    versions           = ch_versions
}
