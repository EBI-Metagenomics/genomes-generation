include { CHECKM } from '../modules/checkm'
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

    main:

    emit:

}




checkm
parse_checkm

echo "genome,completeness,contamination" > all.stats.clean
cut -f1-3 ../checkm_results.tab | sed 's/\s/,/g' | sed 's/\,/\.fa\,/' >> all.stats.clean