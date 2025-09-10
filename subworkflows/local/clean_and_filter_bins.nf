include { CAT as MAG_CLEANUP_CAT    } from '../../modules/local/cat/cat/main'
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

    main:

    ch_versions = Channel.empty()

    bins = input_data.transpose()

    MAG_CLEANUP_CAT ( 
        bins, 
        file(params.cat_db_folder, checkIfExists: true),
        file(params.cat_taxonomy_db, checkIfExists: true),
        file(params.cat_diamond_db, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( MAG_CLEANUP_CAT.out.versions )

    // input: tuple(meta, bin, summary, names)
    // output: tuple(meta, clean.fa)
    DETECT_CONTAMINATION ( 
        MAG_CLEANUP_CAT.out.cat_results 
    )
    ch_versions = ch_versions.mix( DETECT_CONTAMINATION.out.versions )

    GUNC ( 
        DETECT_CONTAMINATION.out.cleaned_fasta, 
        file(params.gunc_db, checkIfExists: true) 
    )
    ch_versions = ch_versions.mix( GUNC.out.versions )

    filtered_bins = GUNC.out.tuple_gunc_result.filter({
        it[2].name.contains('_complete.txt')
    }).map({ _name, cluster_fasta, _cluster_gunc ->
        return cluster_fasta
    })

    emit:
    bins        = filtered_bins
    gunc_report = GUNC.out.gunc_result.collectFile(name: "gunc_report.txt")
    versions    = ch_versions
}
