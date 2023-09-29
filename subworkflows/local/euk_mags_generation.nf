/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGN                          } from './alignment'
include { ALIGN as ALIGN_BINS            } from './alignment'

include { BUSCO                          } from '../../modules/local/busco/main'
include { EUKCC as EUKCC_CONCOCT         } from '../../modules/local/eukcc/main'
include { EUKCC as EUKCC_METABAT         } from '../../modules/local/eukcc/main'
include { LINKTABLE as LINKTABLE_CONCOCT } from '../../modules/local/eukcc/main'
include { LINKTABLE as LINKTABLE_METABAT } from '../../modules/local/eukcc/main'
include { DREP                           } from '../../modules/local/drep/main'
include { DREP as DREP_MAGS              } from '../../modules/local/drep/main'
include { BREADTH_DEPTH                  } from '../../modules/local/breadth_depth/main'
include { BUSCO_EUKCC_QC                 } from '../../modules/local/qc/main'
include { BAT as EUK_TAXONOMY            } from '../../modules/local/CAT/BAT/main'
include { BAT_TAXONOMY_WRITER            } from '../../modules/local/bat_taxonomy_writer/main'


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
    tuple val(meta), path(input_files)
    val output_name

    output:
    path("${output_name}"), emit: modified_result

    script:
    """
    echo "genome,completeness,contamination" > ${output_name}
    cat ${input_files} | grep -v "completeness" | cut -f1-3 | tr '\t' ',' >> ${output_name}
    """
}

process FILTER_QS50 {
    tag "${meta.id}"

    label 'process_light'

    input:
    tuple val(meta), path(quality_file), path(concoct_bins), path(metabat_bins), path(concoct_bins_merged), path(metabat_bins_merged)

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
        for i in \$(ls ${concoct_bins} | grep -w -f filtered_genomes.txt); do
            cp ${concoct_bins}/\${i} output_genomes; done
        for i in \$(ls ${metabat_bins} | grep -w -f filtered_genomes.txt); do
            cp ${metabat_bins}/\${i} output_genomes; done
        for i in \$(ls ${concoct_bins_merged} | grep -w -f filtered_genomes.txt); do
            cp ${concoct_bins_merged}/\${i} output_genomes; done
        for i in \$(ls ${metabat_bins_merged} | grep -w -f filtered_genomes.txt); do
            cp ${metabat_bins_merged}/\${i} output_genomes; done

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
    cat_db
    cat_taxonomy_db

    main:

    ch_versions = Channel.empty()

    /* split the inputs */
    align_input = assemblies_reads_bins.map(item -> tuple(item[0], item[1], item[2]))
    reads = assemblies_reads_bins.map(item -> tuple(item[0], item[2]))
    bins_concoct = assemblies_reads_bins.map(item -> tuple(item[0], item[3]))
    bins_metabat = assemblies_reads_bins.map(item -> tuple(item[0], item[4]))

    ALIGN( align_input )

    ch_versions.mix( ALIGN.out.versions.first() )

    // -- concoct
    binner1 = channel.value("concoct")

    LINKTABLE_CONCOCT( ALIGN.out.assembly_bam.combine(bins_concoct, by: [0]), binner1 ) // output: tuple(meta, links.csv, bin_dir)

    EUKCC_CONCOCT( binner1, LINKTABLE_CONCOCT.out.links_table, eukcc_db )

    // ch_versions.mix( LINKTABLE_CONCOCT.out.versions.first() )
    // ch_versions.mix( EUKCC_CONCOCT.out.versions.first() )

    // -- metabat2
    binner2 = channel.value("metabat2")

    LINKTABLE_METABAT( ALIGN.out.assembly_bam.combine(bins_metabat, by: [0]), binner2 )

    EUKCC_METABAT( binner2, LINKTABLE_METABAT.out.links_table, eukcc_db )

    // ch_versions.mix( LINKTABLE_METABAT.out.versions.first() )
    // ch_versions.mix( EUKCC_METABAT.out.versions.first() )

    // -- prepare quality file
    combine_quality = EUKCC_CONCOCT.out.eukcc_csv.combine( EUKCC_METABAT.out.eukcc_csv, by: [0] )

    // "genome,completeness,contamination" //
    function_cat_csv = { item ->
        def meta = item[0]
        def list_files = [item[1], item[2]]
        return tuple(meta, list_files)
    }

    CONCATENATE_QUALITY_FILES( combine_quality.map(function_cat_csv), channel.value("quality_eukcc.csv") )

    // ch_versions.mix( CONCATENATE_QUALITY_FILES.out.versions.first() )

    quality = CONCATENATE_QUALITY_FILES.out.concatenated_result

    // -- qs50 -- //
    collect_data = quality.combine(bins_concoct, by: 0) \
        .combine(bins_metabat, by: [0]) \
        .combine(EUKCC_CONCOCT.out.eukcc_results, by: [0]) \
        .combine(EUKCC_METABAT.out.eukcc_results, by: [0])

    FILTER_QS50(collect_data)

    // ch_versions.mix( FILTER_QS50.out.versions.first() )

    // input: tuple (meta, genomes/*, quality_file)
    euk_drep_args = channel.value('-pa 0.80 -sa 0.99 -nc 0.40 -cm larger -comp 49 -con 21')

    DREP( FILTER_QS50.out.qs50_filtered_genomes, euk_drep_args, channel.value('euk') )

    // ch_versions.mix( DREP.out.versions.first() )

    // -- aggregate by samples
    quality_all_csv = quality.map{it->it[1]}.collectFile(name: "all.csv", newLine: false)

    MODIFY_QUALITY_FILE( quality_all_csv, channel.value("aggregated_euk_quality.csv"))

    // ch_versions.mix( MODIFY_QUALITY_FILE.out.versions.first() )

    // MODIFY_QUALITY_FILE.out.modified_result.view()

    aggregated_quality = MODIFY_QUALITY_FILE.out.modified_result.map{ it ->
                                                                        def meta = [:]
                                                                        meta.id = "aggregated"
                                                                        return tuple(meta, it)
                                                                    }
    // -- drep MAGs --//
    euk_drep_args_mags = channel.value('-pa 0.80 -sa 0.95 -nc 0.40 -cm larger -comp 49 -con 21')

    combine_drep = DREP.out.dereplicated_genomes.map{item -> item[1] } \
        .flatten() \
        .collect() \
        .map{ it ->
                def meta = [:]
                meta.id = "aggregated"
                return tuple(meta, it)
            }

    DREP_MAGS( combine_drep.combine(aggregated_quality, by: [0]), euk_drep_args_mags, channel.value('euk_mags') )

    // ch_versions.mix( DREP_MAGS.out.versions.first() )

    // -- coverage -- //
    bins_alignment = DREP.out.dereplicated_genomes.combine(reads, by: [0]) // tuple(meta, drep_genomes, [reads]),...

    spreadBins = { record ->
        if (record[1] instanceof List) {
            return record
        } else {
            return tuple(record[0], [record[1]], record[2])
        }
    }

    bins_alignment_by_bins = bins_alignment.map(spreadBins).transpose(by: [1])  // tuple(meta, MAG1, [reads]); tuple(meta, MAG2, [reads])
    
    // TODO: leave only DREP_MAGS instead of DREP
    ALIGN_BINS( bins_alignment_by_bins ) // out: [meta, fasta, bam, bai]

    // ch_versions.mix( ALIGN_BINS.out.versions.first() )

    // input: tuple(meta, MAG, bam, bai)
    BREADTH_DEPTH( ALIGN_BINS.out.assembly_bam )

    // ch_versions.mix( BREADTH_DEPTH.out.versions.first() )

    // -- busco MAGs -- //
    drep_result = DREP_MAGS.out.dereplicated_genomes.map { it -> it[1] }.flatten()

    BUSCO( drep_result, busco_db )

    // ch_versions.mix( BUSCO.out.versions.first() )

    // -- QC MAGs -- //
    BUSCO_EUKCC_QC( aggregated_quality, BUSCO.out.busco_summary.collect(), DREP_MAGS.out.dereplicated_genomes_list.map{ it -> it[1] })

    // ch_versions.mix( BUSCO.out.versions.first() )

    // -- BAT --//
    EUK_TAXONOMY( drep_result, cat_db, cat_taxonomy_db )

    BAT_TAXONOMY_WRITER(EUK_TAXONOMY.out.bat_names.collect())

    // ch_versions.mix( EUK_TAXONOMY.out.versions.first() )
    // ch_versions.mix( BAT_TAXONOMY_WRITER.out.versions.first() )

    emit:
    euk_quality = FILTER_QS50.out.qs50_filtered_genomes.map(item -> tuple(item[0], item[2]))
    drep_output = DREP.out.dereplicated_genomes
    versions = ch_versions
}

