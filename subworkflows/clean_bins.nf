/*
    ~~~~~~~~~~~~~~~~~~
     Collect bins
    ~~~~~~~~~~~~~~~~~~
*/
process COLLECT_BINS {
    container 'quay.io/microbiome-informatics/metawrap:latest'

    input:
    path bins

    output:
    path "Clean_bins/*", emit: chosen_bins

    script:
    """
    mkdir Clean_bins
    cp ${bins}/*fa Clean_bins/
    cd Clean_bins
    CLEAN=$(ls *clean.fa)
    for C in $CLEAN; do mv $C ${C//_clean}; done
    """
}

include { MAG_CLEANUP_CAT } from '../modules/CAT'
include { DETECT_DECONTAMINATION } from '../modules/detect_decontamination'
include { SELECT_SEQS } from '../modules/select_seqs'
/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow CLEAN_BINS {
    take:
        name
        bins
        cat_db
        cat_diamond_db
        cat_taxonomy_db
    main:
        MAG_CLEANUP_CAT(name, cat_db, cat_taxonomy_db, cat_diamond_db, bins)
        DETECT_DECONTAMINATION(name, MAG_CLEANUP_CAT.out.cat_summary, MAG_CLEANUP_CAT.out.cat_names)
        SELECT_SEQS(name, bins, DETECT_DECONTAMINATION.out.cont_contigs)
        //COLLECT_BINS(SELECT_SEQS.out.clean_bins)
    emit:
        cleaned_bins = SELECT_SEQS.out.clean_bins //COLLECT_BINS.out.chosen_bins
}
