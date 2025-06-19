/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EUKCC_MERGE as EUKCC_MERGE_CONCOCT          } from '../subworkflows/local/eukcc_merge'
include { EUKCC_MERGE as EUKCC_MERGE_METABAT          } from '../subworkflows/local/eukcc_merge'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BAT                                         } from '../modules/local/cat/bat/bat'
include { BAT_TAXONOMY_WRITER                         } from '../modules/local/cat/bat/bat_taxonomy_writer'
include { BUSCO                                       } from '../modules/local/busco/busco'
include { BUSCO_EUKCC_QC                              } from '../modules/local/busco/busco_qc'
include { COVERAGE_RECYCLER as COVERAGE_RECYCLER_EUK  } from '../modules/local/coverage_recycler'
include { DREP                                        } from '../modules/local/drep'
include { DREP as DREP_MAGS                           } from '../modules/local/drep'
include { PIGZ as COMPRESS_MAGS                       } from '../modules/local/compress/pigz'
include { PIGZ as COMPRESS_BINS                       } from '../modules/local/compress/pigz'
include { FILTER_QUALITY                              } from '../modules/local/euk_utils'
include { CONCATENATE_QUALITY_FILES                   } from '../modules/local/euk_utils'
include { MODIFY_QUALITY_FILE                         } from '../modules/local/euk_utils'

/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes workflow
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

    /* -- qs50 -- //
    // [meta, concoct_bins, metabat_bins, merged_concoct, merged_metabat]
    // combine concoct, metabat bins with merged bins (if any)
    */
    eukcc_bins = input.concoct_input.map{meta, _assemblies, _reads, concoct, _depth -> [meta, concoct]}
        .join( input.metabat_input.map{meta, _assemblies, _reads, metabat, _depth -> [meta, metabat]} )
        .join( EUKCC_MERGE_CONCOCT.out.eukcc_merged_bins )
        .join( EUKCC_MERGE_METABAT.out.eukcc_merged_bins )
        .map { meta, concoct, metabat, eukcc_concoct, eukcc_metabat ->
            // Collect all files from all folders
            def all_files = []
            [concoct, metabat, eukcc_concoct, eukcc_metabat].each { folder ->
                if (folder && folder.exists()) {
                    if (folder.isDirectory()) {
                        folder.listFiles()?.each { file ->
                            if (file.isFile()) all_files.add(file)
                        }
                    } else {
                        all_files.add(folder)
                    }
                }
            }
            [meta, all_files]
        }

    eukcc_bins.view()

    FILTER_QUALITY( 
        quality.join( eukcc_bins )
    )
    bins = FILTER_QUALITY.out.qs50_filtered_genomes.map{ _meta, genomes, _quality -> genomes }.flatten()

    // input: tuple (meta, genomes/*, quality_file)
    DREP( 
        FILTER_QUALITY.out.qs50_filtered_genomes, 
        params.euk_drep_args, 
        "eukaryotes" 
    )
    ch_versions = ch_versions.mix(DREP.out.versions)

    // -- aggregate by samples
    MODIFY_QUALITY_FILE( 
        quality.map { meta, quality_file -> quality_file }.collectFile(name: "all.csv", newLine: false), 
        "aggregated_euk_quality.csv"
    )
    aggregated_quality = MODIFY_QUALITY_FILE.out.modified_result.map { modified_csv ->
        return tuple([id: "aggregated"], modified_csv)
    }

    // -- drep MAGs --//

    DREP_MAGS( 
        DREP.out.dereplicated_genomes.map{ _meta, drep_genomes -> drep_genomes }.flatten().collect()
          .map{ agg_genomes ->
              return tuple([id: "aggregated"], agg_genomes)
          }.join( aggregated_quality ), 
        params.euk_drep_args_mags, 
        "eukaryotes"
    )
    ch_versions = ch_versions.mix( DREP_MAGS.out.versions)

    // ---- coverage generation ----- //

    COVERAGE_RECYCLER_EUK(
        DREP_MAGS.out.dereplicated_genomes,
        input.concoct_input.map{_meta, _assemblies, _reads, _metabat, depth -> depth}.collectFile(name: "euks_depth.txt.gz")
    )
    ch_versions = ch_versions.mix( COVERAGE_RECYCLER_EUK.out.versions)

    // ---- QC generation----- //
    // -- BUSCO MAG --//
    drep_result = DREP_MAGS.out.dereplicated_genomes.map { _meta, drep_genomes -> drep_genomes }.flatten()

    BUSCO( 
        drep_result, 
        file(params.busco_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( BUSCO.out.versions)

    BUSCO_EUKCC_QC( 
        aggregated_quality.map { _meta, agg_quality_file -> agg_quality_file },
        BUSCO.out.busco_summary.collect(), 
        DREP_MAGS.out.dereplicated_genomes_list.map { _meta, drep_genomes -> drep_genomes }
    )
    ch_versions = ch_versions.mix( BUSCO_EUKCC_QC.out.versions)

    // ---- Taxonomy generation ----- //
    // -- BAT --//
    BAT( 
        drep_result, 
        file(params.cat_db_folder, checkIfExists: true), 
        file(params.cat_taxonomy_db, checkIfExists: true) 
    )
    ch_versions = ch_versions.mix( BAT.out.versions)

    BAT_TAXONOMY_WRITER( 
        BAT.out.bat_names.collect() 
    )
    ch_versions = ch_versions.mix( BAT_TAXONOMY_WRITER.out.versions)

    // compress euk genomes
    COMPRESS_MAGS(
        drep_result
    )
    // publish genomes
    COMPRESS_MAGS.out.compressed.subscribe({ cluster_fasta ->
        cluster_fasta.copyTo("${params.outdir}/genomes_drep/eukaryotes/genomes/${cluster_fasta.name}")
    })

    // compress euk bins
    COMPRESS_BINS(
        bins
    )
    // publish bins
    COMPRESS_BINS.out.compressed.subscribe({ cluster_fasta ->
        cluster_fasta.copyTo("${params.outdir}/bins/eukaryotes/${cluster_fasta.name.split('_')[0]}/${cluster_fasta.name}")
    })

    // logging collection
    ch_log = ch_log.mix (EUKCC_MERGE_METABAT.out.progress_log )
    ch_log = ch_log.mix (EUKCC_MERGE_CONCOCT.out.progress_log )
    ch_log = ch_log.mix( FILTER_QUALITY.out.progress_log )
    ch_log = ch_log.mix( DREP.out.progress_log )
    ch_log = ch_log.mix( DREP_MAGS.out.progress_log )

    /* 
    // ---- Generate files for bins upload ----
    */
    if ( params.upload_bins ) {
        // choose bins that are not MAGs
        drep_result.view()
        bins.view()
        // COVERAGE_RECYCLER_EUK
        // BAT
        // BAT_TAXONOMY_WRITER
        // PREPARE_UPLOAD_FILES
    }


    emit:
    genomes                   = COMPRESS_MAGS.out.compressed.collect()
    stats                     = BUSCO_EUKCC_QC.out.eukcc_final_qc
    coverage                  = COVERAGE_RECYCLER_EUK.out.mag_coverage.map{ _meta, coverage_file -> coverage_file }.collect()
    taxonomy                  = BAT_TAXONOMY_WRITER.out.all_bin2classification
    versions                  = ch_versions
    progress_log              = ch_log
    // for multiqc
    samtools_idxstats_metabat = EUKCC_MERGE_METABAT.out.samtools_idxstats
    samtools_idxstats_concoct = EUKCC_MERGE_CONCOCT.out.samtools_idxstats
    busco_short_summary       = BUSCO.out.busco_summary
}