/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EUKCC_MERGE as EUKCC_MERGE_CONCOCT          } from '../subworkflows/local/eukcc_merge'
include { EUKCC_MERGE as EUKCC_MERGE_METABAT          } from '../subworkflows/local/eukcc_merge'

/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow EUK_MAGS_GENERATION {
    take:
    assemblies_and_reads
    concoct_bins
    metabat_bins
    jgi_depth

    main:

    ch_versions = Channel.empty()
    ch_log      = Channel.empty()

    /*
    * LINKTABLE and EUKCC for euk bin refinement for CONCOCT bins
    * input: tuple( meta, assembly_file, [raw_reads], bins, depths )
    */
    EUKCC_MERGE_CONCOCT(
        assemblies_and_reads.join(concoct_bins).join(jgi_depth)
    )
    ch_versions = ch_versions.mix( EUKCC_MERGE_CONCOCT.out.versions )

    /*
    * LINKTABLE and EUKCC for euk bin refinement for METABAT bins
    * input: tuple( meta, assembly_file, [raw_reads], bins, depths )
    */
    EUKCC_MERGE_METABAT(
        assemblies_and_reads.join(metabat_bins).join(jgi_depth)
    )
    ch_versions = ch_versions.mix( EUKCC_MERGE_METABAT.out.versions )

    emit:
    versions                  = ch_versions
}