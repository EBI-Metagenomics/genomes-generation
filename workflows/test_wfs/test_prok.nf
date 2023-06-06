/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/
ref_genome = file(params.ref_genome, checkIfExists: true)
ref_genome_name = file(params.ref_genome_name)
ref_cat_diamond = file("${params.CAT_ref_db}/${params.cat_diamond_db_name}", checkIfExists: true)
ref_catdb = file("${params.CAT_ref_db}/${params.cat_db_name}", checkIfExists: true)
ref_cat_taxonomy = file("${params.CAT_ref_db}/${params.cat_taxonomy_db}", checkIfExists: true)
ref_eukcc = file("${params.eukcc_ref_db}", checkIfExists: true)
ref_gunc = file("${params.gunc_ref_db}", checkIfExists: true)
ref_checkm = file("${params.checkm_ref_db}", checkIfExists: true)
ref_rfam_rrna_models = file("${params.rfam_rrna_models}", checkIfExists: true)
ref_gtdbtk = file("${params.gtdbtk}", checkIfExists: true)

/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/
bin_ref_1 = channel.fromPath("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover2/data/bin_ref_SRR6311991/*")
bin_ref_2 = channel.fromPath("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover2/data/bin_ref_SRR6312008/*")
/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { CLEAN_AND_FILTER_BINS } from '../subworkflows/clean_and_filter_bins'
include { CHECKM2 } from '../modules/checkm2'
include { DREP } from '../modules/drep'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {

    data = bin_ref_1.concat(bin_ref_2)

    groupReads = { fastq ->
        def cluster = fastq.toString().tokenize("/")[-2].tokenize("/")[0]
        return tuple(cluster, fastq)
    }

    input_data = data.map(groupReads).groupTuple()
    CLEAN_AND_FILTER_BINS(input_data, ref_catdb.first(), ref_cat_diamond.first(), ref_cat_taxonomy.first(), ref_gunc.first())

    CHECKM2(CLEAN_AND_FILTER_BINS.out.filtered_bins, ref_checkm.first())

    drep_input = CLEAN_AND_FILTER_BINS.out.filtered_bins.combine(CHECKM2.out.checkm_table, by: 0)
    prok_drep_args = channel.value('-pa 0.9 -sa 0.95 -nc 0.6 -cm larger -comp 50 -con 5')
    DREP(drep_input, prok_drep_args)
}
