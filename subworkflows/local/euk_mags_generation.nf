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
    cat ${quality_table_csv} | grep -v "completeness" | cut -f1-3 | tr '\t' ',' >> ${output_name}
    """
}

process FILTER_QS50 {
    tag "${meta.id}"

    label 'process_light'

    input:
    tuple val(meta), path(quality_file), path(concoct_bins, stageAs: "concoct_bins/*"), path(metabat_bins, stageAs: "metabat_bins/*"), path(concoct_bins_merged, stageAs: "concoct_bins_merged/*"), path(metabat_bins_merged, stageAs: "metabat_bins_merged/*")

    output:
    tuple val(meta), path("output_genomes/*"), path("quality_file.csv"), emit: qs50_filtered_genomes, optional: true

    script:
    """
    # prepare drep quality file
    cat ${quality_file} | grep -v "completeness" |\
        awk '{{if(\$2 - 5*\$3 >=50){{print \$0}}}}' |\
        sort -k 2,3 -n | cut -f1 > filtered_genomes.txt
    BINS=\$(cat filtered_genomes.txt | wc -l)
    mkdir -p output_genomes
    if [ \$BINS -lt 2 ];
    then
        touch quality_file.csv
    else
        for i in \$(ls concoct_bins | grep -w -f filtered_genomes.txt); do
            cp concoct_bins/\${i} output_genomes; done
        for i in \$(ls metabat_bins | grep -w -f filtered_genomes.txt); do
            cp metabat_bins/\${i} output_genomes; done
        for i in \$(ls concoct_bins_merged | grep -w -f filtered_genomes.txt); do
            cp concoct_bins_merged/\${i} output_genomes; done
        for i in \$(ls metabat_bins_merged | grep -w -f filtered_genomes.txt); do
            cp metabat_bins_merged/\${i} output_genomes; done

        echo "genome,completeness,contamination" > quality_file.csv
        grep -w -f filtered_genomes.txt ${quality_file} | cut -f1-3 | tr '\\t' ',' >> quality_file.csv
    fi
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

    /* split the inputs */
    assemblies_reads_bins.multiMap { meta, assembly, reads, concoct_folder, metabat_folder -> 
        assembly_and_reads: [ meta, assembly, reads ]
        reads: [ meta, reads ]
        bins_concoct: [ meta, concoct_folder ]
        bins_metabat: [ meta, metabat_folder ]
    }.set {
        input
    }

    ALIGN( input.assembly_and_reads )

    ch_versions.mix( ALIGN.out.versions.first() )

    // -- concoct -- //
    binner1 = channel.value("concoct")

    LINKTABLE_CONCOCT( ALIGN.out.assembly_bam.join( input.bins_concoct ), binner1 ) // output: tuple(meta, links.csv)

    EUKCC_CONCOCT( binner1, LINKTABLE_CONCOCT.out.links_table.join ( input.bins_concoct ), eukcc_db )

    ch_versions.mix( LINKTABLE_CONCOCT.out.versions.first() )
    ch_versions.mix( EUKCC_CONCOCT.out.versions.first() )

    // -- metabat2
    binner2 = channel.value("metabat2")

    LINKTABLE_METABAT( ALIGN.out.assembly_bam.join( input.bins_metabat ), binner2 )

    metabat_linktable_bins = LINKTABLE_METABAT.out.links_table.join( input.bins_metabat ).filter { meta, link, bins -> bins.size() > 0 }

    EUKCC_METABAT( binner2, metabat_linktable_bins, eukcc_db )

    ch_versions.mix( LINKTABLE_METABAT.out.versions.first() )
    ch_versions.mix( EUKCC_METABAT.out.versions.first() )

    // -- prepare quality file
    combine_quality = EUKCC_CONCOCT.out.eukcc_csv.join( EUKCC_METABAT.out.eukcc_csv )

    // "genome,completeness,contamination" //
    functionCATCSV = { item ->
        def meta = item[0]
        def list_files = [item[1], item[2]]
        return tuple(meta, list_files)
    }

    CONCATENATE_QUALITY_FILES( combine_quality.map( functionCATCSV ), channel.value("quality_eukcc.csv") )

    quality = CONCATENATE_QUALITY_FILES.out.concatenated_result

    // -- qs50 -- //
    collect_data = quality.join( input.bins_concoct ) \
        .join( input.bins_metabat ) \
        .join( EUKCC_CONCOCT.out.eukcc_results ) \
        .join( EUKCC_METABAT.out.eukcc_results )

    FILTER_QS50( collect_data )

    // input: tuple (meta, genomes/*, quality_file)
    DREP( FILTER_QS50.out.qs50_filtered_genomes, params.euk_drep_args, channel.value('eukaryotes') )

    ch_versions.mix( DREP.out.versions.first() )

    // -- aggregate by samples
    // TODO: this collectFile is incorrect
    quality_all_csv = quality.map { meta, quality_file -> quality_file }.collectFile(name: "all.csv", newLine: false)

    MODIFY_QUALITY_FILE( quality_all_csv, channel.value("aggregated_euk_quality.csv"))

    // ch_versions.mix( MODIFY_QUALITY_FILE.out.versions.first() )

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

    DREP_MAGS( combine_drep.join( aggregated_quality ), params.euk_drep_args_mags, channel.value('aggregated_eukaryotes') )

    ch_versions.mix( DREP_MAGS.out.versions.first() )

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

    ch_versions.mix( ALIGN_BINS.out.versions.first() )

    // ---- coverage generation ----- //
    ch_summarizedepth_input = ALIGN_BINS.out.assembly_bam.map { meta, assembly, bams, bais ->
        [ meta, bams, bais ]
    }
    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS ( ch_summarizedepth_input )
    COVERAGE_RECYCLER_EUK(
        DREP_MAGS.out.dereplicated_genomes,
        METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth.map{ meta, depth -> depth }.collect()
    )

    ch_versions.mix( METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions.first() )
    ch_versions.mix( COVERAGE_RECYCLER_EUK.out.versions.first() )

    // ---- QC generation----- //
    // -- BUSCO MAG --//
    drep_result = DREP_MAGS.out.dereplicated_genomes.map { meta, drep_genomes -> drep_genomes }.flatten()

    BUSCO( drep_result, busco_db )

    ch_versions.mix( BUSCO.out.versions.first() )

    BUSCO_EUKCC_QC( 
        aggregated_quality.map { meta, agg_quality_file -> agg_quality_file },
        BUSCO.out.busco_summary.collect(), 
        DREP_MAGS.out.dereplicated_genomes_list.map { meta, drep_genomes -> drep_genomes }
    )

    ch_versions.mix( BUSCO_EUKCC_QC.out.versions.first() )

    // ---- Taxonomy generation ----- //
    // -- BAT --//
    BAT( drep_result, cat_db_folder, cat_taxonomy_db )

    BAT_TAXONOMY_WRITER( BAT.out.bat_names.collect() )

    ch_versions.mix( BAT.out.versions.first() )
    ch_versions.mix( BAT_TAXONOMY_WRITER.out.versions.first() )

    emit:
    euk_quality = BUSCO_EUKCC_QC.out.busco_final_qc
    taxonomy = BAT_TAXONOMY_WRITER.out.all_bin2classification
    genomes = drep_result
    versions = ch_versions
}
