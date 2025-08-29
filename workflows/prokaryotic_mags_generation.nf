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

include { BINETTE                    } from '../modules/ebi-metagenomics/binette/main'
include { CHECKM2                    } from '../modules/local/checkm2/main'
include { DREP_DEREPLICATE           } from '../modules/nf-core/drep/dereplicate/main'
include { PIGZ as COMPRESS_MAGS      } from '../modules/local/compress/pigz'
include { PIGZ as COMPRESS_BINS      } from '../modules/local/compress/pigz'
include { COVERAGE_RECYCLER          } from '../modules/local/coverage_recycler/main'
include { DETECT_RRNA                } from '../modules/local/detect_rrna/main'
include { GTDBTK                     } from '../modules/local/gtdbtk/main'
include { GTDBTK_TO_NCBI_TAXONOMY    } from '../modules/local/gtdbtk/gtdb_to_ncbi_majority_vote/main'
include { PROPAGATE_TAXONOMY_TO_BINS } from '../modules/local/propagate_taxonomy_to_bins/main'

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
    all_bins = CLEAN_AND_FILTER_BINS.out.bins.collect().map { bins_list ->
        [ [id: "aggregated" ], bins_list ]
    }

    /* --  Estimate completeness and contamination for all project bins -- */
    CHECKM2 (
        all_bins,
        file(params.checkm2_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( CHECKM2.out.versions.first() )

    // TODO add FILTER_QUALITY step here before dereplication

    /* --  Dereplicate and filter good quality bins -- */
    DREP_DEREPLICATE (
        CHECKM2.out.bins_and_stats,
        [[id:''], []]   // No previous dRep work directory
    )
    ch_versions = ch_versions.mix( DREP_DEREPLICATE.out.versions )

    /* --  Calculate coverage for all genomes -- */
    COVERAGE_RECYCLER ( 
        all_bins,
        all_metabat_depths
    )
    ch_versions = ch_versions.mix( COVERAGE_RECYCLER.out.versions.first() )

    /* --  Detect RNA for cluster representatives -- */
    dereplicated_genomes = DREP_DEREPLICATE.out.fastas.map { _meta, bins_list -> bins_list }.flatten()
    DETECT_RRNA(
        dereplicated_genomes,
        file(params.rfam_rrna_models, checkIfExists: true)
    )
    rna_out = Channel.empty()
    rna_out = rna_out.mix( DETECT_RRNA.out.rrna_out_files.collect(), DETECT_RRNA.out.trna_out_files.collect()  )
    ch_versions = ch_versions.mix( DETECT_RRNA.out.versions.first() )

    /* --  Assign taxonomy to cluster representatives -- */
    GTDBTK (
        dereplicated_genomes.collect(),
        file(params.gtdbtk_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( GTDBTK.out.versions.first() )

    /* --  Convert GTDB taxonomy to NCBI -- */
    GTDBTK_TO_NCBI_TAXONOMY (
        GTDBTK.out.gtdbtk_output,
        file(params.gtdbtk_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( GTDBTK_TO_NCBI_TAXONOMY.out.versions.first() )

    /* --  Propagate taxonomy to from cluster representatives to cluster members -- */
    // TODO what will happen if all clusters are singletons
    clustering_csvs = DREP_DEREPLICATE.out.summary_tables
        .map { _meta, summary_table -> 
            summary_table.findAll { table -> 
                !(table.name in ["Cdb.csv", "Wdb.csv"] )
            } 
        }
    PROPAGATE_TAXONOMY_TO_BINS (
        DREP_DEREPLICATE.out.summary_tables,
        clustering_csvs
    )
    ch_versions = ch_versions.mix( PROPAGATE_TAXONOMY_TO_BINS.out.versions.first() )

    /* --  Compress MAGs -- */
    COMPRESS_MAGS (
        dereplicated_genomes
    )
    compressed_mags = COMPRESS_MAGS.out.compressed.subscribe({ cluster_fasta ->
        cluster_fasta.copyTo("${params.outdir}/${params.subdir_proks}/${params.subdir_mags}/${cluster_fasta.name}")
    })

    // TODO bins

    /* --  Compress bins -- */
    COMPRESS_BINS (
        CLEAN_AND_FILTER_BINS.out.bins
    )

    /* --  Finalize logging -- */
    ch_log = Channel.empty()
    ch_log = ch_log.mix( BINETTE.out.progress_log )
    ch_log = ch_log.mix( CHECKM2.out.progress_log )
    ch_log = ch_log.mix( DREP_DEREPLICATE.out.progress_log )

    emit:
    mags_fastas  = COMPRESS_MAGS.out.compressed.collect()
    bins_fastas  = COMPRESS_BINS.out.compressed.collect()
    stats        = CHECKM2.out.bins_and_stats.map { _map, _bins, stats -> stats }
    coverage     = COVERAGE_RECYCLER.out.mag_coverage.map{ _meta, coverage_file -> coverage_file }.collect()
    mags_rna     = rna_out.collect()
    taxonomy     = PROPAGATE_TAXONOMY_TO_BINS.out.ncbi_taxonomy.map { _meta, taxonomy_file -> taxonomy_file }.collect()
    versions     = ch_versions
    progress_log = ch_log
}