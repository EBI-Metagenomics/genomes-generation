/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGN                          } from './local/subwf_alignment'
include { ALIGN as ALIGN_BINS            } from './local/subwf_alignment'
include { BUSCO                          } from '../modules/local/busco/main'
include { EUKCC as EUKCC_CONCOCT         } from '../modules/local/eukcc/main'
include { EUKCC as EUKCC_METABAT         } from '../modules/local/eukcc/main'
include { EUKCC_MAG                      } from '../modules/local/eukcc/main'
include { EUKCC_SINGLE                   } from '../modules/local/eukcc/main'
include { LINKTABLE as LINKTABLE_CONCOCT } from '../modules/local/eukcc/main'
include { LINKTABLE as LINKTABLE_METABAT } from '../modules/local/eukcc/main'
include { DREP                           } from '../modules/local/drep/main'
include { DREP_MAGS                      } from '../modules/local/drep/main'
include { BREADTH_DEPTH                  } from '../modules/local/breadth_depth/main'
include { BUSCO_QC                       } from '../modules/local/busco_qc/main'
include { BAT as EUK_TAXONOMY            } from '../modules/local/CAT/BAT/main'
include { BAT_TAXONOMY_WRITER            } from '../modules/local/bat_taxonomy_writer/main'


process CONCATENATE_QUALITY_FILES {
    tag "${meta.id}"

    input:
    tuple val(meta), path(input_files)
    val output_name

    output:
    tuple val(meta), path("${output_name}"), emit: concatenated_result

    script:
    """
    echo "bin\tcompleteness\tcontamination\tncbi_lng" > ${output_name}
    cat ${input_files} | grep -v "completeness" >> ${output_name}
    """
}

process MODIFY_QUALITY_FILE {

    input:
    path(input_file)
    val output_name

    output:
    path("${output_name}"), emit: modified_result

    script:
    """
    echo "genome,completeness,contamination" > ${output_name}
    cat ${input_file} | grep -v "completeness" | cut -f1-3 | tr '\t' ',' >> ${output_name}
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


workflow EUK_SUBWF {
    take:
        input_data  // tuple( meta, assembly_file, [raw_reads], concoct_folder, metabat_folder )
        eukcc_db
        busco_db
        cat_db
        cat_taxonomy_db

    main:
        align_input = input_data.map(item -> tuple(item[0], item[1], item[2]))
        reads = input_data.map(item -> tuple(item[0], item[2]))
        bins_concoct = input_data.map(item -> tuple(item[0], item[3]))
        bins_metabat = input_data.map(item -> tuple(item[0], item[4]))

        ALIGN(align_input)

        // -- concoct
        binner1 = channel.value("concoct")
        LINKTABLE_CONCOCT(ALIGN.out.output.combine(bins_concoct, by: 0), binner1)       // output: tuple(meta, links.csv, bin_dir)
        EUKCC_CONCOCT(binner1, LINKTABLE_CONCOCT.out.links_table, eukcc_db.first())

        // -- metabat2
        binner2 = channel.value("metabat2")
        LINKTABLE_METABAT(ALIGN.out.output.combine(bins_metabat, by: 0), binner2)
        EUKCC_METABAT(binner2, LINKTABLE_METABAT.out.links_table, eukcc_db.first())

        // -- prepare quality file
        combine_quality = EUKCC_CONCOCT.out.eukcc_csv.combine(EUKCC_METABAT.out.eukcc_csv, by: 0)
        // "genome,completeness,contamination"
        functionCatCSV = { item ->
            def meta = item[0]
            def list_files = [item[1], item[2]]
            //def combine_quality_file = list_files.collectFile(name: "quality_eukcc.csv", newLine: true)
            return tuple(meta, list_files)
        }
        CONCATENATE_QUALITY_FILES(combine_quality.map(functionCatCSV), channel.value("quality_eukcc.csv"))
        quality = CONCATENATE_QUALITY_FILES.out.concatenated_result

        // -- qs50
        collect_data = quality.combine(bins_concoct, by: 0).combine(bins_metabat, by: 0).combine(EUKCC_CONCOCT.out.eukcc_results, by: 0).combine(EUKCC_METABAT.out.eukcc_results, by: 0)
        FILTER_QS50(collect_data)

        size_filtered_bins = FILTER_QS50(collect_data).out.qs50_filtered_genomes.map{it->it[1]}.collect().size()
        size_refined_bins.subscribe { println "Euk QS50: $it.value" }
        // -- drep
        // input: tuple (meta, genomes/*, quality_file)
        euk_drep_args = channel.value('-pa 0.80 -sa 0.99 -nc 0.40 -cm larger -comp 49 -con 21')
        DREP(FILTER_QS50.out.qs50_filtered_genomes, euk_drep_args, channel.value('euk'))

        // -- coverage
        bins_alignment = DREP.out.dereplicated_genomes.combine(reads, by:0)  // tuple(meta, drep_genomes, [reads]),...
        spreadBins = { record ->
            if (record[1] instanceof List) {
                return record}
            else {
                return tuple(record[0], [record[1]], record[2])}
        }
        bins_alignment_by_bins = bins_alignment.map(spreadBins).transpose(by:1)  // tuple(meta, MAG1, [reads]); tuple(meta, MAG2, [reads])
        ALIGN_BINS(bins_alignment_by_bins)      // out: tuple(meta, [bam, bai])
        breadth_input = ALIGN_BINS.out.output.combine(bins_alignment_by_bins, by:0).map(item -> tuple(item[0], item[1], item[2]))
        // input: tuple(meta, [bam, bai], MAG)
        BREADTH_DEPTH(breadth_input)

        // -- aggregate by samples
        combine_drep = DREP.out.dereplicated_genomes.map(item -> item[1]).flatten().collect()
        agg_quality = quality.map(item -> item[1]).flatten().collectFile(name:"aggregated_quality.txt", skip:1, keepHeader:true)
        MODIFY_QUALITY_FILE(agg_quality, channel.value("euk_quality.csv"))
        // -- drep MAGs
        euk_drep_args_mags = channel.value('-pa 0.80 -sa 0.95 -nc 0.40 -cm larger -comp 49 -con 21')
        DREP_MAGS(channel.value("aggregated"), combine_drep, MODIFY_QUALITY_FILE.out.modified_result, euk_drep_args_mags, channel.value('euk_mags'))

        // -- eukcc MAGs
        //EUKCC_SINGLE(DREP_MAGS.out.dereplicated_genomes.flatten(), eukcc_db)

        // -- busco MAGs - Varsha
        BUSCO(DREP_MAGS.out.dereplicated_genomes.flatten(), busco_db.first())

        // -- QC MAGs - Ales
	BUSCO_QC(CONCATENATE_QUALITY_FILES.out.concatenated_result.map{it[1]}, BUSCO.out.busco_summary.collect())

	// -- BAT - Ales
	EUK_TAXONOMY(DREP_MAGS.out.dereplicated_genomes.flatten(), cat_db.first(), cat_taxonomy_db.first())
	BAT_TAXONOMY_WRITER(EUK_TAXONOMY.out.bat_names.collect())


    emit:
        euk_quality = FILTER_QS50.out.qs50_filtered_genomes.map(item -> tuple(item[0], item[2]))
        drep_output = DREP.out.dereplicated_genomes

}

