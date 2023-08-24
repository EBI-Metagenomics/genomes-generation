include { CAT as MAG_CLEANUP_CAT    } from '../../modules/local/CAT/CAT/main'
include { DETECT_CONTAMINATION      } from '../../modules/local/detect_contamination/main'
include { GUNC                      } from '../../modules/local/gunc/main'
/*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Run subworkflow cleaning and filtering with GUNC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run subwf on each bin.fa
    input: bins - channel of one bin.fa
*/
workflow CLEAN_AND_FILTER_BINS {
    take:
        input_data // meta, [bin.fa, ...]
        cat_db
        cat_diamond_db
        cat_taxonomy_db
        gunc_db
    main:
        bins = input_data.transpose()
        MAG_CLEANUP_CAT(bins, cat_db, cat_taxonomy_db, cat_diamond_db)

        // input: tuple(meta, bin, summary, names)
        // output: tuple(meta, clean.fa)
        DETECT_CONTAMINATION(MAG_CLEANUP_CAT.out.cat_results)

        GUNC(DETECT_CONTAMINATION.out.cleaned_fasta, gunc_db.first())
        filtered_bins = GUNC.out.tuple_gunc_result.filter({
                it[2].name.contains('_complete.txt')
            }).map({ name, cluster_fasta, cluster_gunc ->
                return cluster_fasta
            })
        filtered_bins.view()
    emit:
        bins = filtered_bins
        gunc_report = GUNC.out.gunc_result.collectFile(name: "gunc_report.txt")
}
