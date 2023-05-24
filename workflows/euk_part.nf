/*
    ~~~~~~~~~~~~~~~~~~~~~~~~
     Eukaryotes subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EUKCC as EUKCC_CONCOCT } from '../modules/eukcc'
include { EUKCC as EUKCC_METABAT } from '../modules/eukcc'
include { LINKTABLE as LINKTABLE_CONCOCT } from '../modules/eukcc'
include { LINKTABLE as LINKTABLE_METABAT } from '../modules/eukcc'
//include { FILTER_QS50 } from '../modules/filter_qs50_euk'

workflow EUK_SUBWF {
    take:
        bins_concoct    // tuple: name, folder
        bins_metabat    // tuple: name, folder
        bam             // tuple: name, [bam, bam.bai]
        eukcc_db
    main:
        // concoct
        binner1 = channel.value("concoct")
        LINKTABLE_CONCOCT(bins_concoct, bam)
        EUKCC_CONCOCT(binner1, LINKTABLE_CONCOCT.out.links_table, eukcc_db.first(), bins_concoct)

        // metabat2
        binner2 = channel.value("metabat2")
        LINKTABLE_METABAT(bins_metabat, bam)
        EUKCC_METABAT(binner2, LINKTABLE_METABAT.out.links_table, eukcc_db.first(), bins_metabat)

        // prepare quality file
        combine_quality = EUKCC_CONCOCT.out.eukcc_csv.combine(EUKCC_METABAT.out.eukcc_csv, by: 0)
        combine_quality.view()
        //quality = EUKCC_CONCOCT.out.eukcc_csv.map{
        //    item -> item[1].concat(item[2]).collectFile(name: "quality_eukcc.csv").text().into {
        //        file -> file.write("genome,completeness,contamination\n") }
        //        }

        //FILTER_QS50(quality, bins_concoct, bins_metabat, EUKCC_CONCOCT.out.eukcc_results, EUKCC_METABAT.out.eukcc_results)
        //DREP()
        // aggregate outputs
        // coverage
        // drep MAGs
        // eukcc MAGs
        // busco MAGs
        // QC MAGs
    emit:
        euk_quality = EUKCC_CONCOCT.out.eukcc_csv
}