include { MAG_CLEANUP_CAT } from '../modules/CAT'
include { DETECT_DECONTAMINATION } from '../modules/detect_decontamination'
include { SELECT_SEQS } from '../modules/select_seqs'
/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
    Run subwf on each bin.fa
    input: bins - channel of bin.fa files
*/
workflow CLEAN_BINS {
    take:
        name
        bins
        cat_db
        cat_diamond_db
        cat_taxonomy_db
    main:
        MAG_CLEANUP_CAT(name, cat_db.first(), cat_taxonomy_db.first(), cat_diamond_db.first(), bins)
        DETECT_DECONTAMINATION(name, MAG_CLEANUP_CAT.out.cat_summary, MAG_CLEANUP_CAT.out.cat_names)
        SELECT_SEQS(name, bins, DETECT_DECONTAMINATION.out.cont_contigs)
        //COLLECT_BINS
    emit:
        cleaned_bins = SELECT_SEQS.out.clean_bins
}
