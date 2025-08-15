/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLEAN_AND_FILTER_BINS              } from '../subworkflows/local/clean_and_filter_bins'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BINETTE                            } from '../modules/ebi-metagenomics/binette/main'
include { CHECKM2                            } from '../modules/local/checkm2/main'
include { DREP                               } from '../modules/local/drep/main'
include { PIGZ as COMPRESS_MAGS              } from '../modules/local/compress/pigz'
include { COVERAGE_RECYCLER                  } from '../modules/local/coverage_recycler/main'
include { DETECT_RRNA                        } from '../modules/local/detect_rrna/main'
include { GTDBTK                             } from '../modules/local/gtdbtk/main'
include { GTDBTK_TO_NCBI_TAXONOMY            } from '../modules/local/gtdbtk/gtdb_to_ncbi_majority_vote/main'

include { CHECKM2_TABLE_FOR_DREP_GENOMES     } from '../modules/local/utils'
include { CHANGE_UNDERSCORE_TO_DOT           } from '../modules/local/utils'

/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Prokaryotes workflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROK_MAGS_GENERATION {

    take:

    collected_binners_assembly_and_depth // tuple( meta, concoct, metabat, maxbin, assembly_fasta, depth_file)

    main:

    ch_versions = Channel.empty()

    all_metabat_depths = collected_binners_assembly_and_depth.map { _meta, _concoct, _metabat, _maxbin, _assembly_fasta, depth_file -> depth_file }.collect()

    /* --  Bins refinement -- */
    ch_binette_input = collected_binners_assembly_and_depth.map { meta, concoct_bins, maxbin_bins, metabat_bins, assembly, _depth ->
        def valid_bins = [concoct_bins, maxbin_bins, metabat_bins].findAll { it != null && it != [] }
        [
            meta,
            valid_bins,    // Combined binning results
            assembly,
            []             // Placeholder for predicted proteins (optional BINETTE input)
        ]
    }

    BINETTE(
        ch_binette_input,
        "fasta",
        file(params.checkm2_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( BINETTE.out.versions )

    /* --  Clean bins -- */
    CLEAN_AND_FILTER_BINS(
        BINETTE.out.refined_bins
    )
    ch_versions = ch_versions.mix( CLEAN_AND_FILTER_BINS.out.versions.first() )

    /* --  aggregate all bins  -- */
    all_bins = CLEAN_AND_FILTER_BINS.out.bins.collect().map { it ->
        def meta = [:]
        meta.id = "aggregated"
        return tuple( meta, it )
    }

    /* --  Estimate completeness and contamination for all project bins -- */
    CHECKM2 (
        "aggregated",
        all_bins,
        file(params.checkm2_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( CHECKM2.out.versions.first() )

    /* --  Dereplicate and filter good quality bins -- */
    prok_drep_args = channel.value('-pa 0.9 -sa 0.95 -nc 0.6 -cm larger -comp 50 -con 5')

    DREP ( 
        CHECKM2.out.stats, 
        prok_drep_args, 
        "prokaryotes" 
    )
    dereplicated_genomes = DREP.out.dereplicated_genomes.map { it -> it[1] }.flatten()
    ch_versions = ch_versions.mix( DREP.out.versions.first() )

    /* --  Calculate coverage -- */
    COVERAGE_RECYCLER ( 
        DREP.out.dereplicated_genomes, 
        all_metabat_depths
    )
    ch_versions = ch_versions.mix( COVERAGE_RECYCLER.out.versions.first() )

    /* --  Change underscores to dots in fasta headers for MAGs -- */
    CHANGE_UNDERSCORE_TO_DOT ( 
        dereplicated_genomes 
    )

    /* --  Detect RNA -- */
    DETECT_RRNA(
        dereplicated_genomes,
        file(params.rfam_rrna_models, checkIfExists: true)
    )
    rna_out = Channel.empty()
    rna_out = rna_out.mix( DETECT_RRNA.out.rrna_out_files.collect(), DETECT_RRNA.out.trna_out_files.collect()  )
    ch_versions = ch_versions.mix( DETECT_RRNA.out.versions.first() )

    /* --  Assign taxonomy -- */
    GTDBTK (
        CHANGE_UNDERSCORE_TO_DOT.out.return_files.collect(),
        file(params.gtdbtk_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( GTDBTK.out.versions.first() )

    /* --  Convert GTDB taxonomy to NCBI -- */
    GTDBTK_TO_NCBI_TAXONOMY (
        GTDBTK.out.gtdbtk_output,
        file(params.gtdbtk_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( GTDBTK_TO_NCBI_TAXONOMY.out.versions.first() )

    /* --  Modify CheckM table -- */
    // checkm_results_mags.txt 
    // Both channels will have only one element
    CHECKM2_TABLE_FOR_DREP_GENOMES (
        CHECKM2.out.stats.map { _map, _bins, stats -> stats },
        DREP.out.dereplicated_genomes_list.map { _meta, genomes_list_tsv -> genomes_list_tsv }
    )

    /* --  Compress MAGs -- */
    COMPRESS_MAGS (
        CHANGE_UNDERSCORE_TO_DOT.out.return_files.collect().flatten()
    )
    compressed_bins = COMPRESS_MAGS.out.compressed.subscribe({ cluster_fasta ->
        cluster_fasta.copyTo("${params.outdir}/genomes_drep/prokaryotes/genomes/${cluster_fasta.name}")
    })

    /* --  Finalize logging -- */
    ch_log = Channel.empty()
    ch_log = ch_log.mix( BINETTE.out.progress_log )
    ch_log = ch_log.mix( CHECKM2.out.progress_log )
    ch_log = ch_log.mix( DREP.out.progress_log )

    emit:

    genomes = COMPRESS_MAGS.out.compressed.collect()
    stats = CHECKM2_TABLE_FOR_DREP_GENOMES.out.checkm_results_mags
    coverage = COVERAGE_RECYCLER.out.mag_coverage.map{ _meta, coverage_file -> coverage_file }.collect()
    rna = rna_out.collect()
    taxonomy = GTDBTK_TO_NCBI_TAXONOMY.out.ncbi_taxonomy
    versions = ch_versions
    progress_log = ch_log
}