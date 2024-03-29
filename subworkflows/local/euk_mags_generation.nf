/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGN                                      } from './alignment'
include { ALIGN as ALIGN_BINS                        } from './alignment'
include { BUSCO                                      } from '../../modules/local/busco/main'
include { EUKCC as EUKCC_CONCOCT                     } from '../../modules/local/eukcc/main'
include { EUKCC as EUKCC_METABAT                     } from '../../modules/local/eukcc/main'
include { LINKTABLE as LINKTABLE_CONCOCT             } from '../../modules/local/eukcc/main'
include { LINKTABLE as LINKTABLE_METABAT             } from '../../modules/local/eukcc/main'
include { DREP                                       } from '../../modules/local/drep/main'
include { DREP as DREP_MAGS                          } from '../../modules/local/drep/main'
include { BUSCO_EUKCC_QC                             } from '../../modules/local/qc/main'
include { BAT                                        } from '../../modules/local/cat/bat/main'
include { BAT_TAXONOMY_WRITER                        } from '../../modules/local/bat_taxonomy_writer/main'
include { COVERAGE_RECYCLER as COVERAGE_RECYCLER_EUK } from '../../modules/local/coverage_recycler/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS       } from '../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { GZIP as GZIP_MAGS                          } from '../../modules/local/utils'
include { GZIP as GZIP_BINS                          } from '../../modules/local/utils'


process CONCATENATE_QUALITY_FILES {
    tag "${meta.id}"

    input:
    tuple val(meta), path(input_files)
    val output_name

    output:
    tuple val(meta), path("${meta.id}.${output_name}"), emit: concatenated_result

    script:
    """
    echo "bin\tcompleteness\tcontamination\tncbi_lng" > "${meta.id}.${output_name}"
    for i in ${input_files}; do
        tail -n +2 \${i} > help_file
        if [ -s help_file ]; then
            grep -v "completeness" \${i} >> "${meta.id}.${output_name}"
        fi
    done
    rm help_file
    """
}

process MODIFY_QUALITY_FILE {

    input:
    path(quality_table_csv)
    val output_name

    output:
    path("${output_name}"), emit: modified_result

    script:
    """
    echo "genome,completeness,contamination" > ${output_name}
    grep -v "completeness" ${quality_table_csv} | cut -f1-3 | tr '\t' ',' >> ${output_name} || true
    """
}

process FILTER_QUALITY {
    tag "${meta.id}"

    label 'process_light'

    input:
    tuple val(meta), path(quality_file), path(concoct_bins, stageAs: "concoct_bins/*"), path(metabat_bins, stageAs: "metabat_bins/*"), path(concoct_bins_merged, stageAs: "concoct_bins_merged/*"), path(metabat_bins_merged, stageAs: "metabat_bins_merged/*")

    output:
    tuple val(meta), path("output_genomes/*"), path("quality_file.csv"), emit: qs50_filtered_genomes, optional: true
    path "progress.log"                                       , emit: progress_log

    script:
    """
    mkdir -p output_genomes
    touch quality_file.csv

    echo "Prepare drep quality"
    # prepare drep quality file
    grep -v "completeness" ${quality_file} |\
    awk '{{if(\$2>=50 && \$2<=100 && \$3>=0 && \$3<=5){{print \$0}}}}' |\
    sort -k 2,3 -n | cut -f1 > filtered_genomes.txt || true

    echo "bins count"
    export BINS=\$(cat filtered_genomes.txt | wc -l)
    echo "\$BINS"
    if [ \$BINS -lt 2 ];
    then
        echo "No genomes"
    else
        for i in \$(ls concoct_bins | grep -w -f filtered_genomes.txt); do
            cp concoct_bins/\${i} output_genomes; done
        for i in \$(ls metabat_bins | grep -w -f filtered_genomes.txt); do
            cp metabat_bins/\${i} output_genomes; done
        for i in \$(ls concoct_bins_merged/merged_bins | grep -w -f filtered_genomes.txt); do
            cp concoct_bins_merged/merged_bins/\${i} output_genomes; done
        for i in \$(ls metabat_bins_merged/merged_bins | grep -w -f filtered_genomes.txt); do
            cp metabat_bins_merged/merged_bins/\${i} output_genomes; done

        echo "genome,completeness,contamination" > quality_file.csv
        grep -w -f filtered_genomes.txt ${quality_file} | cut -f1-3 | tr '\\t' ',' >> quality_file.csv
    fi

    cat <<-END_LOGGING > progress.log
    ${meta.id}\t${task.process}
        concoct_bins: \$(ls concoct_bins | wc -l), metabat_bins: \$(ls metabat_bins | wc -l), concoct_bins_merged: \$(ls concoct_bins_merged/merged_bins | wc -l), metabat_bins_merged: \$(ls metabat_bins_merged/merged_bins | wc -l), output_genomes: \$(ls output_genomes | wc -l)
    END_LOGGING
    """
}


workflow EUK_MAGS_GENERATION {
    take:
    assemblies_reads_bins  // tuple( meta, assembly_file, [raw_reads], concoct_folder, metabat_folder )
    eukcc_db
    busco_db
    cat_db_folder
    cat_diamond_db
    cat_taxonomy_db

    main:

    ch_versions = Channel.empty()
    ch_log      = Channel.empty()

    /* split the inputs */
    assemblies_reads_bins.multiMap { meta, assembly, reads, concoct_bins, metabat_bins ->
        assembly_and_reads: [ meta, assembly, reads ]
        reads: [ meta, reads ]
        bins_concoct: [ meta, concoct_bins ]
        bins_metabat: [ meta, metabat_bins ]
    }.set {
        input
    }

    ALIGN( input.assembly_and_reads )

    ch_versions = ch_versions.mix( ALIGN.out.versions )

    // -- concoct -- //
    binner1 = channel.value("concoct")

    LINKTABLE_CONCOCT( ALIGN.out.assembly_bam.join( input.bins_concoct ), binner1 ) // output: tuple(meta, links.csv)

    EUKCC_CONCOCT( binner1, LINKTABLE_CONCOCT.out.links_table.join ( input.bins_concoct ), eukcc_db )

    ch_versions = ch_versions.mix( LINKTABLE_CONCOCT.out.versions.first() )
    ch_versions = ch_versions.mix( EUKCC_CONCOCT.out.versions.first() )

    // -- metabat2
    binner2 = channel.value("metabat2")

    LINKTABLE_METABAT( ALIGN.out.assembly_bam.join( input.bins_metabat ), binner2 )

    metabat_linktable_bins = LINKTABLE_METABAT.out.links_table.join( input.bins_metabat ).filter { meta, link, bins -> bins.size() > 0 }

    EUKCC_METABAT( binner2, metabat_linktable_bins, eukcc_db )

    ch_versions = ch_versions.mix( LINKTABLE_METABAT.out.versions )
    ch_versions = ch_versions.mix( EUKCC_METABAT.out.versions )

    // -- prepare quality file
    combine_quality = EUKCC_CONCOCT.out.eukcc_csv.join( EUKCC_METABAT.out.eukcc_csv )

    // "genome,completeness,contamination" //
    functionCATCSV = { item ->
        def meta = item[0]
        def list_files = [item[1], item[2]]
        return tuple(meta, list_files)
    }

    CONCATENATE_QUALITY_FILES( combine_quality.map( functionCATCSV ), "quality_eukcc.csv" )

    quality = CONCATENATE_QUALITY_FILES.out.concatenated_result

    // -- qs50 -- //
    // [meta, concoct_bins, metabat_bins, merged_concoct, merged_metabat]
    // combine concoct, metabat bins with merged bins (if any)
    collect_data = quality.join( input.bins_concoct ) \
        .join( input.bins_metabat ) \
        .join( EUKCC_CONCOCT.out.eukcc_merged_bins ) \
        .join( EUKCC_METABAT.out.eukcc_merged_bins )

    FILTER_QUALITY( collect_data )

    // input: tuple (meta, genomes/*, quality_file)
    DREP( FILTER_QUALITY.out.qs50_filtered_genomes, params.euk_drep_args, "eukaryotes" )

    ch_versions = ch_versions.mix( DREP.out.versions.first() )

    // -- aggregate by samples
    quality_all_csv = quality.map { meta, quality_file -> quality_file }.collectFile(name: "all.csv", newLine: false)

    MODIFY_QUALITY_FILE( quality_all_csv, "aggregated_euk_quality.csv")

    aggregated_quality = MODIFY_QUALITY_FILE.out.modified_result.map { modified_csv ->
        return tuple([id: "aggregated"], modified_csv)
    }

    // -- drep MAGs --//

    combine_drep = DREP.out.dereplicated_genomes.map{ meta, drep_genomes -> drep_genomes } \
        .flatten() \
        .collect() \
        .map{ agg_genomes ->
            return tuple([id: "aggregated"], agg_genomes)
        }

    DREP_MAGS( combine_drep.join( aggregated_quality ), params.euk_drep_args_mags, 'eukaryotes' )

    ch_versions = ch_versions.mix( DREP_MAGS.out.versions.first() )

    // -- coverage -- //
    bins_alignment = DREP.out.dereplicated_genomes.join( input.reads ) // tuple(meta, drep_genomes, [reads]),...

    spreadBins = { record ->
        if (record[1] instanceof List) {
            return record
        } else {
            return tuple(record[0], [record[1]], record[2])
        }
    }

    bins_alignment_by_bins = bins_alignment.map( spreadBins ).transpose(by: [1])  // tuple(meta, MAG1, [reads]); tuple(meta, MAG2, [reads])

    ALIGN_BINS( bins_alignment_by_bins ) // out: [meta, fasta, bam, bai]

    ch_versions = ch_versions.mix( ALIGN_BINS.out.versions )

    // ---- coverage generation ----- //
    ch_summarizedepth_input = ALIGN_BINS.out.assembly_bam.map { meta, assembly, bams, bais ->
        [ meta, bams, bais ]
    }
    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )
    euks_depth = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth.map{ meta, depth -> depth }.collectFile(name: "euks_depth.txt.gz")
    COVERAGE_RECYCLER_EUK(
        DREP_MAGS.out.dereplicated_genomes,
        euks_depth
    )

    ch_versions = ch_versions.mix( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first() )
    ch_versions = ch_versions.mix( COVERAGE_RECYCLER_EUK.out.versions.first() )

    // ---- QC generation----- //
    // -- BUSCO MAG --//
    drep_result = DREP_MAGS.out.dereplicated_genomes.map { meta, drep_genomes -> drep_genomes }.flatten()

    BUSCO( drep_result, busco_db )

    ch_versions = ch_versions.mix( BUSCO.out.versions.first() )

    BUSCO_EUKCC_QC( 
        aggregated_quality.map { meta, agg_quality_file -> agg_quality_file },
        BUSCO.out.busco_summary.collect(), 
        DREP_MAGS.out.dereplicated_genomes_list.map { meta, drep_genomes -> drep_genomes }
    )

    ch_versions = ch_versions.mix( BUSCO_EUKCC_QC.out.versions.first() )

    // ---- Taxonomy generation ----- //
    // -- BAT --//
    BAT( drep_result, cat_db_folder, cat_taxonomy_db )

    BAT_TAXONOMY_WRITER( BAT.out.bat_names.collect() )

    ch_versions = ch_versions.mix( BAT.out.versions.first() )
    ch_versions = ch_versions.mix( BAT_TAXONOMY_WRITER.out.versions.first() )

    // compress euk genomes
    GZIP_MAGS(drep_result)
    compressed_genomes = GZIP_MAGS.out.compressed.subscribe({ cluster_fasta ->
        cluster_fasta.copyTo("${params.outdir}/genomes_drep/eukaryotes/genomes/${cluster_fasta.name}")
    })
    // compress euk bins
    bins_to_compress = FILTER_QUALITY.out.qs50_filtered_genomes.map{ meta, genomes, quality -> genomes }
    GZIP_BINS(bins_to_compress.flatten())
    compressed_bins = GZIP_BINS.out.compressed.subscribe({ cluster_fasta ->
        cluster_fasta.copyTo("${params.outdir}/bins/eukaryotes/${cluster_fasta.name.split('_')[0]}/${cluster_fasta.name}")
    })

    ch_log = ch_log.mix (EUKCC_METABAT.out.progress_log )
    ch_log = ch_log.mix (EUKCC_CONCOCT.out.progress_log )
    ch_log = ch_log.mix( FILTER_QUALITY.out.progress_log )
    ch_log = ch_log.mix( DREP.out.progress_log )
    ch_log = ch_log.mix( DREP_MAGS.out.progress_log )

    emit:
    genomes = GZIP_MAGS.out.compressed.collect()
    stats = BUSCO_EUKCC_QC.out.eukcc_final_qc
    coverage = COVERAGE_RECYCLER_EUK.out.mag_coverage.map{ meta, coverage_file -> coverage_file }.collect()
    taxonomy = BAT_TAXONOMY_WRITER.out.all_bin2classification
    versions = ch_versions
    progress_log = ch_log
}
