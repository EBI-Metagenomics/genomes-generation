include { CHECKM2 } from '../modules/checkm2'
include { PARSE_CHECKM } from '../modules/parse_checkm'
/*
    ~~~~~~~~~~~~~~~~~~~~~~
     Run checkm subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CHECKM_SUBWF {
    take:
        name
        bins
        checkm_db
    main:
        CHECKM2(bins, checkm_db)
        //PARSE_CHECKM(CHECKM2.out.bins_taxonomy, CHECKM2.out.bins_qa)
    emit:
        files_checkm = CHECKM2.out.checkm_output
}

//echo "genome,completeness,contamination" > all.stats.clean
//cut -f1-3 ../checkm_results.tab | sed 's/\s/,/g' | sed 's/\,/\.fa\,/' >> all.stats.clean