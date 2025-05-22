/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGNMENT_LINKTABLE as LINKTABLE           } from '../../modules/local/bwa_mem2/align_linktable/main'
include { EUKCC as EUKCC                             } from '../../modules/local/eukcc/main'

/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     EukCC subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EUKCC_MERGE {
    take:
    assemblies_reads_bins_depth  // tuple( meta, assembly_file, [raw_reads], bins, depths )

    main:

    ch_versions = Channel.empty()
    ch_log      = Channel.empty()

    /* split the inputs */
    assemblies_reads_bins_depth.multiMap { meta, assembly, reads, bin_folder, depths ->
        assembly: [ meta, assembly ]
        reads: [ meta, reads ]
        bin_folder: [ meta, bin_folder ]
        depths: [ meta, depths ]
    }.set {
        input
    }

    /*
    * eukcc linktable.py
    * output: tuple(meta, links.csv)
    */
    LINKTABLE(
        input.reads.join(input.assembly).join(input.bin_folder).join(input.depths)
    )
    ch_versions = ch_versions.mix( LINKTABLE.out.versions )


    EUKCC (
        LINKTABLE.out.links_table.join( input.bin_folder ),
        file(params.eukcc_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( EUKCC.out.versions)
    ch_log = ch_log.mix (EUKCC.out.progress_log )


    emit:
    eukcc_csv                 = EUKCC.out.eukcc_csv
    eukcc_merged_bins         = EUKCC.out.eukcc_merged_bins
    samtools_idxstats         = LINKTABLE.out.idxstats
    versions                  = ch_versions
    progress_log              = ch_log
}