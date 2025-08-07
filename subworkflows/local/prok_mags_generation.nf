/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Prokaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLEAN_AND_FILTER_BINS              } from './clean_and_filter_bins'

include { BINETTE                            } from '../../modules/ebi-metagenomics/binette/main'
include { CHECKM2                            } from '../../modules/local/checkm2/main'
include { DREP                               } from '../../modules/local/drep/main'
include { COVERAGE_RECYCLER                  } from '../../modules/local/coverage_recycler/main'
include { DETECT_RRNA                        } from '../../modules/local/detect_rrna/main'
include { GTDBTK                             } from '../../modules/local/gtdbtk/main'
include { GTDBTK_TO_NCBI_TAXONOMY            } from '../../modules/local/gtdbtk/gtdb_to_ncbi_majority_vote/main'
include { CHANGE_UNDERSCORE_TO_DOT           } from '../../modules/local/utils'
include { GZIP as GZIP_MAGS                  } from '../../modules/local/utils'

process CHECKM2_TABLE_FOR_DREP_GENOMES {

    input:
    path(checkm_filtered_genomes_dir)
    path(dereplicated_genomes_tsv)

    output:
    path("checkm_results_mags.tab"), emit: checkm_results_mags

    script:
    """
    grep -f ${dereplicated_genomes_tsv} ${checkm_filtered_genomes_dir} > checkm_results_mags.tab || true
    """
}


workflow PROK_MAGS_GENERATION {

    take:
    collected_binners_assembly_and_depth // tuple( meta, concoct, metabat, maxbin, assembly_fasta, depth_file)
    cat_db_folder
    cat_diamond_db
    cat_taxonomy_db
    gunc_db
    checkm2_db
    gtdbtk_db
    rfam_rrna_models

    main:

    ch_versions = Channel.empty()
    ch_log = Channel.empty()

    ch_binette_input = collected_binners_assembly_and_depth.map { meta, concoct_bins, maxbin_bins, metabat_bins, assembly, _ -> 
        def valid_bins = [concoct_bins, maxbin_bins, metabat_bins].findAll { it != null && it != [] }
        [
            meta,
            valid_bins,    // Combined binning results
            assembly,
            []             // Placeholder for predicted proteins (optional BINETTE input)
        ]
    }

    all_metabat_depths = collected_binners_assembly_and_depth.map { meta, _concoct, _metabat, _maxbin, _assembly_fasta, depth_file -> depth_file }.collect()

    // -- bin refinement //
    BINETTE( ch_binette_input, "fasta", checkm2_db )
    ch_versions = ch_versions.mix( BINETTE.out.versions )

    // -- clean bins
    CLEAN_AND_FILTER_BINS( 
        BINETTE.out.refined_bins,
        cat_db_folder,
        cat_diamond_db,
        cat_taxonomy_db,
        gunc_db
    )

    ch_versions = ch_versions.mix( CLEAN_AND_FILTER_BINS.out.versions.first() )

    // -- aggregate bins by samples
    // -- checkm2 on ALL bins in all samples
    all_bins = CLEAN_AND_FILTER_BINS.out.bins.collect().map { it ->
        def meta = [:]
        meta.id = "aggregated"
        return tuple( meta, it )
    }

    CHECKM2( "aggregated", all_bins, checkm2_db )

    ch_versions = ch_versions.mix( CHECKM2.out.versions.first() )

    // -- drep
    prok_drep_args = channel.value('-pa 0.9 -sa 0.95 -nc 0.6 -cm larger -comp 50 -con 5')

    DREP( CHECKM2.out.stats, prok_drep_args, "prokaryotes" )

    ch_versions = ch_versions.mix( DREP.out.versions.first() )

    // -- coverage -- //
    COVERAGE_RECYCLER( DREP.out.dereplicated_genomes, all_metabat_depths)

    ch_versions = ch_versions.mix( COVERAGE_RECYCLER.out.versions.first() )

    dereplicated_genomes = DREP.out.dereplicated_genomes.map { it -> it[1] }.flatten()

    CHANGE_UNDERSCORE_TO_DOT( dereplicated_genomes )

    // -- RNA -- //
    DETECT_RRNA( 
        dereplicated_genomes,
        rfam_rrna_models
    )
    rna_out = Channel.empty()
    rna_out = rna_out.mix( DETECT_RRNA.out.rrna_out_files.collect(), DETECT_RRNA.out.trna_out_files.collect()  )
    ch_versions = ch_versions.mix( DETECT_RRNA.out.versions.first() )

    // -- Taxonomy --//
    GTDBTK( CHANGE_UNDERSCORE_TO_DOT.out.return_files.collect(), gtdbtk_db )
    ch_versions = ch_versions.mix( GTDBTK.out.versions.first() )
    GTDBTK_TO_NCBI_TAXONOMY(GTDBTK.out.gtdbtk_output, gtdbtk_db)
    ch_versions = ch_versions.mix( GTDBTK_TO_NCBI_TAXONOMY.out.versions.first() )

    // -- checkm_results_mags.txt -- //
    // Both channels will have only one element
    CHECKM2_TABLE_FOR_DREP_GENOMES(
        CHECKM2.out.stats.map { map, bins, stats -> stats },
        DREP.out.dereplicated_genomes_list.map { meta, genomes_list_tsv -> genomes_list_tsv }
    )
    // compress prok genomes
    GZIP_MAGS(CHANGE_UNDERSCORE_TO_DOT.out.return_files.collect().flatten())
    compressed_bins = GZIP_MAGS.out.compressed.subscribe({ cluster_fasta ->
        cluster_fasta.copyTo("${params.outdir}/genomes_drep/prokaryotes/genomes/${cluster_fasta.name}")
    })

    ch_log = ch_log.mix( BINETTE.out.progress_log )
    ch_log = ch_log.mix( CHECKM2.out.progress_log )
    ch_log = ch_log.mix( DREP.out.progress_log )

    emit:
    genomes = GZIP_MAGS.out.compressed.collect()
    stats = CHECKM2_TABLE_FOR_DREP_GENOMES.out.checkm_results_mags
    coverage = COVERAGE_RECYCLER.out.mag_coverage.map{ meta, coverage_file -> coverage_file }.collect()
    rna = rna_out.collect()
    taxonomy = GTDBTK_TO_NCBI_TAXONOMY.out.ncbi_taxonomy
    versions = ch_versions
    progress_log = ch_log
}
