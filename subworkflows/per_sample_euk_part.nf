/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGN } from '../subworkflows/subwf_alignment'
include { ALIGN as ALIGN_BINS} from '../subworkflows/subwf_alignment'
include { EUKCC as EUKCC_CONCOCT } from '../modules/eukcc'
include { EUKCC as EUKCC_METABAT } from '../modules/eukcc'
include { LINKTABLE as LINKTABLE_CONCOCT } from '../modules/eukcc'
include { LINKTABLE as LINKTABLE_METABAT } from '../modules/eukcc'
include { DREP } from '../modules/drep'
include { DREP as DREP_MAGS } from '../modules/drep'
include { BREADTH_DEPTH } from '../modules/breadth_depth'


process CONCATENATE_QUALITY_FILES {
    tag "${name}"

    publishDir(
        path: "${params.outdir}/qs50",
        mode: "copy"
    )

    input:
    tuple val(name), path(input_files)
    val output_name

    output:
    tuple val(name), path("${output_name}"), emit: concatenated_result

    script:
    """
    echo "bin\tcompleteness\tcontamination\tncbi_lng" > ${output_name}
    cat ${input_files} | grep -v "completeness" >> ${output_name}
    """
}

process FILTER_QS50 {
    tag "${name}"

    publishDir(
        path: "${params.outdir}/qs50",
        mode: "copy"
    )

    label 'process_light'

    input:
    tuple val(name), path(quality_file), path(concoct_bins), path(metabat_bins), path(concoct_bins_merged), path(metabat_bins_merged)

    output:
    tuple val(name), path("output_genomes/*"), path("quality_file.csv"), emit: qs50_filtered_genomes, optional: true

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
        input_data  // tuple( run_accession, assembly_file, [raw_reads], concoct_folder, metabat_folder )
        eukcc_db
    main:

        align_input = input_data.map(item -> tuple(item[0], item[1], item[2]))
        reads = input_data.map(item -> tuple(item[0], item[2]))
        bins_concoct = input_data.map(item -> tuple(item[0], item[3]))
        bins_metabat = input_data.map(item -> tuple(item[0], item[4]))

        ALIGN(align_input)

        // concoct
        binner1 = channel.value("concoct")
        LINKTABLE_CONCOCT(ALIGN.out.annotated_bams.combine(bins_concoct, by: 0), binner1)       // output: tuple(name, links.csv, bin_dir)
        EUKCC_CONCOCT(binner1, LINKTABLE_CONCOCT.out.links_table, eukcc_db.first())

        // metabat2
        binner2 = channel.value("metabat2")
        LINKTABLE_METABAT(ALIGN.out.annotated_bams.combine(bins_metabat, by: 0), binner2)
        EUKCC_METABAT(binner2, LINKTABLE_METABAT.out.links_table, eukcc_db.first())

        // prepare quality file
        combine_quality = EUKCC_CONCOCT.out.eukcc_csv.combine(EUKCC_METABAT.out.eukcc_csv, by: 0)
        // "genome,completeness,contamination"
        functionCatCSV = { item ->
            def name = item[0]
            def list_files = [item[1], item[2]]
            //def combine_quality_file = list_files.collectFile(name: "quality_eukcc.csv", newLine: true)
            return tuple(name, list_files)
        }
        CONCATENATE_QUALITY_FILES(combine_quality.map(functionCatCSV), channel.value("quality_eukcc.csv"))
        quality = CONCATENATE_QUALITY_FILES.out.concatenated_result

        collect_data = quality.combine(bins_concoct, by: 0).combine(bins_metabat, by: 0).combine(EUKCC_CONCOCT.out.eukcc_results, by: 0).combine(EUKCC_METABAT.out.eukcc_results, by: 0)
        FILTER_QS50(collect_data)

        // input: tuple (name, genomes/*, quality_file)
        euk_drep_args = channel.value('-pa 0.80 -sa 0.99 -nc 0.40 -cm larger -comp 49 -con 21')
        DREP(FILTER_QS50.out.qs50_filtered_genomes, euk_drep_args, channel.value('euk'))

        // coverage
        bins_alignment = DREP.out.dereplicated_genomes.combine(reads, by:0).transpose(by:1)  // tuple(name, [bins], [reads])
        ALIGN_BINS(bins_alignment)
        BREADTH_DEPTH(ALIGN_BINS.out.annotated_bams)

        // aggregate outputs
        combine_drep = DREP.out.dereplicated_genomes.map(item -> tuple(channel.value("aggregated"), item[1])).groupTuple()
        combine_drep.view()
        euk_drep_args_mags = channel.value('-pa 0.80 -sa 0.95 -nc 0.40 -cm larger -comp 49 -con 21')
        //DREP_MAGS(combine_drep, euk_drep_args_mags, channel.value('euk_mags'))

        // drep MAGs
        // eukcc MAGs
        // busco MAGs
        // QC MAGs
    emit:
        euk_quality = FILTER_QS50.out.qs50_filtered_genomes.map(item -> tuple(item[0], item[2]))
        drep_output = DREP.out.dereplicated_genomes
}