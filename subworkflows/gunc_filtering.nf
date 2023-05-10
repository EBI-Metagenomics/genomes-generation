include { GUNC } from '../modules/gunc'
/*
    ~~~~~~~~~~~~~~~~~~~~~~
     Run GUNC subworkflow
    ~~~~~~~~~~~~~~~~~~~~~~
*/
workflow FILTER_BINS {
    take:
        name
        bins
    main:
        GUNC()
    emit:
        filtered_bins = GUNC.out.cluster_gunc_result.filter({
                it[2].name.contains('_complete.txt')
            }).map({ cluster_name, cluster_fasta, cluster_gunc ->
                return tuple(cluster_name, cluster_fasta)
            })
}
