/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/
sample_name = channel.value(params.sample_name)
mode = channel.value(params.mode)
contigs = channel.fromPath(params.contigs, checkIfExists: true)
bam = channel.fromPath(params.bam, checkIfExists: true)
bam_index = channel.fromPath(params.bam_index, checkIfExists: true)

bins = channel.fromPath("${params.bins}/*.fa", checkIfExists: true)

if ( params.mode == "paired" ) {
    chosen_reads = channel.fromFilePairs(["${params.paired_end_forward}", "${params.paired_end_reverse}"], checkIfExists: true).map { it[1] }
}
else if ( params.mode == "single" ) {
    chosen_reads = channel.fromPath("${params.single_end}", checkIfExists: true)
}

/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/
ref_genome = channel.fromPath(params.ref_genome, checkIfExists: true)
ref_genome_name = channel.value(params.ref_genome_name)
ref_cat_diamond = channel.fromPath("${params.CAT_ref_db}/${params.cat_diamond_db_name}", checkIfExists: true)
ref_catdb = channel.fromPath("${params.CAT_ref_db}/${params.cat_db_name}", checkIfExists: true)
ref_cat_taxonomy = channel.fromPath("${params.CAT_ref_db}/${params.cat_taxonomy_db}", checkIfExists: true)
ref_gunc = channel.fromPath("${params.gunc_ref_db}", checkIfExists: true)
/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/

include { CLEAN_AND_FILTER_BINS } from '../subworkflows/clean_and_filter_bins'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {

    CLEAN_AND_FILTER_BINS(sample_name, bins, ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc)

}
