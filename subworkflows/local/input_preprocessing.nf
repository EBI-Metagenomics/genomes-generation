//
// Change reads and assemblies input files if necessary
//

include { CHECK_AND_MODIFY_READS } from '../../modules/local/check_and_modify_reads'
include { MODIFY_CONTIGS         } from '../../modules/local/modify_contigs'

include { FASTQC                 } from '../../modules/nf-core/fastqc/main'


workflow INPUT_PREPROCESSING {

    take:
    input_data  // tuple( meta, assembly_file, [raw_reads] )

    main:

    ch_versions = Channel.empty()

    reads = input_data.map { meta, _assembly, reads -> [meta, reads] }
    contigs = input_data.map { meta, assembly, _reads_item -> [meta, assembly] }

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
        * --- MODIFY READS with SEQKIT ( change run_accession to assembly_accession and apply sanity check )
        *  input: tuple(meta, [reads]])
        *  output: tuple(meta, [modified_pe_reads]/modified_se_reads])
        */
        CHECK_AND_MODIFY_READS(
            reads
        )
        ch_versions = ch_versions.mix( CHECK_AND_MODIFY_READS.out.versions.first() )

        // this process returns list of reads if they are PE and path for SE
        // it is required for next step FASTP that SE reads should be also in list
        reads_changed = CHECK_AND_MODIFY_READS.out.modified_reads.map{ meta, reads_item ->
                                                                result = [meta]
                                                                if (!meta.single_end) {
                                                                    result.add(reads_item) }
                                                                else {
                                                                    result.add([reads_item]) }
                                                                return result
                                                              }

        result = MODIFY_CONTIGS.out.underscore_contigs.join(reads_changed)
    }

    emit:
    assembly_and_reads = result // tuple( meta, assembly_file, [raw_reads] )
    fastqc_gzip        = FASTQC.out.zip
    versions           = ch_versions
}
