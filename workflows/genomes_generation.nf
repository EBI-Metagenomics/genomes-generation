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
ref_eukcc = channel.fromPath("${params.eukcc_ref_db}", checkIfExists: true)
ref_gunc = channel.fromPath("${params.gunc_ref_db}", checkIfExists: true)
ref_checkm = channel.fromPath("${params.checkm_ref_db}", checkIfExists: true)
ref_rfam_rrna_models = channel.fromPath("${params.rfam_rrna_models}", checkIfExists: true)
ref_gtdbtk = channel.fromPath("${params.gtdbtk}", checkIfExists: true)
/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_INPUT } from '../subworkflows/prepare_input_files'
include { BINNING } from '../subworkflows/binning'
include { EUK_SUBWF } from '../subworkflows/euk_part'
include { PROK_SUBWF } from '../subworkflows/prok_part'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {

    PREPARE_INPUT(mode, sample_name, contigs, chosen_reads, ref_genome, ref_genome_name)

    BINNING(mode, sample_name, PREPARE_INPUT.out.contigs_fixed, PREPARE_INPUT.out.reads_cleaned)

    PROK_SUBWF(sample_name,
        ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc, ref_checkm, ref_gtdbtk, ref_rfam_rrna_models)

    EUK_SUBWF(sample_name)
}
