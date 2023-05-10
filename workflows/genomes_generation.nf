/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/
sample_name = channel.value(params.sample_name)
mode = channel.value(params.mode)
contigs = channel.fromPath(params.contigs, checkIfExists: true)

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
/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_INPUT } from '../subworkflows/prepare_input_files'
include { BINNING } from '../subworkflows/binning'
include { CLEAN_BINS } from '../subworkflows/clean_bins'
include { FILTER_BINS } from '../subworkflows/gunc_filtering'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {

    PREPARE_INPUT(mode, sample_name, contigs, chosen_reads, ref_genome, ref_genome_name)

    BINNING(mode, sample_name, PREPARE_INPUT.out.contigs_fixed, PREPARE_INPUT.out.reads_cleaned)

    //CLEAN_BINS
    //FILTER_BINS
    //CHECKM_SUBWF
    //EUKCC_SUBWF


}
