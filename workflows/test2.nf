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

bins = channel.fromPath(params.bins, checkIfExists: true)
metabat_depth = channel.fromPath(params.metabat_depth, checkIfExists: true)

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
include { CHECKM2 } from '../modules/checkm2'
include { DREP } from '../modules/drep'
include { DETECT_RRNA } from '../modules/detect_rrna'
include { RENAME } from '../modules/rename'
include { COVERAGE_RECYCLER } from '../modules/cov_recycler'
include { CHANGE_UNDERSCORE_TO_DOT } from '../modules/prepare_data'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {

    CHECKM2(bins, ref_checkm)

    nc_drep = channel.value(0.60)
    DREP(bins, CHECKM2.out.checkm_table, nc_drep)

    RENAME(DREP.out.dereplicated_genomes)

    //DETECT_RRNA(DREP.out.dereplicated_genomes, ref_rfam_rrna_models.first())

    COVERAGE_RECYCLER(DREP.out.dereplicated_genomes, metabat_depth)

    //CHANGE_UNDERSCORE_TO_DOT(DREP.out.dereplicated_genomes)
    //GTDBTK(DREP.out.dereplicated_genomes, ref_gtdbtk)
}
