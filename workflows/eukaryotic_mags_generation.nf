/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EUKCC_MERGE as EUKCC_MERGE_CONCOCT          } from '../subworkflows/local/eukcc_merge'
include { EUKCC_MERGE as EUKCC_MERGE_METABAT          } from '../subworkflows/local/eukcc_merge'
include { CONCATENATE_QUALITY_FILES                   } from '../modules/local/euk_utils'
/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow EUK_MAGS_GENERATION {
    take:
    input_data

    main:

    ch_versions = Channel.empty()
    ch_log      = Channel.empty()

    input_data
        .multiMap { meta, assemblies, reads, concoct, metabat, depth -> 
            concoct_input: [meta, assemblies, reads, concoct, depth]
            metabat_input: [meta, assemblies, reads, metabat, depth]
        }
        .set { input }

    /*
    * LINKTABLE and EUKCC for euk bin refinement for CONCOCT bins
    * input: tuple( meta, assembly_file, [raw_reads], bins, depths )
    */
    EUKCC_MERGE_CONCOCT(
        input.concoct_input
    )
    ch_versions = ch_versions.mix( EUKCC_MERGE_CONCOCT.out.versions )

    /*
    * LINKTABLE and EUKCC for euk bin refinement for METABAT bins
    * input: tuple( meta, assembly_file, [raw_reads], bins, depths )
    */
    EUKCC_MERGE_METABAT(
        input.metabat_input
    )
    ch_versions = ch_versions.mix( EUKCC_MERGE_METABAT.out.versions )

    // "genome,completeness,contamination" //
    functionCATCSV = { item ->
        def meta = item[0]
        def list_files = [item[1], item[2]]
        return tuple(meta, list_files)
    }

    CONCATENATE_QUALITY_FILES(
        EUKCC_MERGE_CONCOCT.out.eukcc_csv
        .join( EUKCC_MERGE_METABAT.out.eukcc_csv )
        .map( functionCATCSV ), "quality_eukcc.csv" 
    )
    quality = CONCATENATE_QUALITY_FILES.out.concatenated_result



    emit:
    versions                  = ch_versions
}