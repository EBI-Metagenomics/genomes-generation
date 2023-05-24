include { MAG_CLEANUP_CAT } from '../modules/CAT'
include { DETECT_DECONTAMINATION } from '../modules/detect_decontamination'
include { SELECT_SEQS } from '../modules/select_seqs'
include { GUNC } from '../modules/gunc'
/*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Run subworkflow cleaning and filtering with GUNC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run subwf on each bin.fa
    input: bins - channel of bin.fa files
*/
workflow CLEAN_AND_FILTER_BINS {
    take:
        bins
        cat_db
        cat_diamond_db
        cat_taxonomy_db
        gunc_db
    main:
        MAG_CLEANUP_CAT(bins, cat_db.first(), cat_taxonomy_db.first(), cat_diamond_db.first())
        DETECT_DECONTAMINATION(MAG_CLEANUP_CAT.out.cat_summary, MAG_CLEANUP_CAT.out.cat_names)
        SELECT_SEQS(bins, DETECT_DECONTAMINATION.out.cont_contigs)
        GUNC(SELECT_SEQS.out.clean_bins, gunc_db.first())
    emit:
        // todo return name
        filtered_bins = GUNC.out.tuple_gunc_result.filter({
                it[1].name.contains('_complete.txt')
            }).map({ cluster_fasta, cluster_gunc ->
                return cluster_fasta
            })
        filtered_bins.view()

        gunc_report = GUNC.out.gunc_result.collectFile(name: "gunc_report.txt")
}
