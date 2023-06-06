/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/
accession = channel.value("SRR6311991")

/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/

/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { BIN_REFINEMENT } from '../../modules/metawrap'
include { CLEAN_AND_FILTER_BINS } from '../../subworkflows/subwf_clean_and_filter_bins'

/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/

collect_binners = channel.from([
                    tuple(
                        "ERZ_1",
                        file("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover/data/bins/concoct_bins"),
                        file("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover/data/bins/maxbin2_bins"),
                        file("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover/data/bins/metabat2_bins")
                        ),
                    tuple(
                        "ERZ_2",
                        file("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover/data/bins/concoct_bins"),
                        file("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover/data/bins/maxbin2_bins"),
                        file("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover/data/bins/metabat2_bins")
                        )
                    ])

metabat_depth_for_coverage = channel.from([
                    tuple(
                        "ERZ_1",
                        file("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover/data/bins/concoct_bins/bin.0.fa")),
                    tuple(
                        "ERZ_2",
                        file("/hps/nobackup/rdf/metagenomics/service-team/users/kates/eukrecover/data/bins/concoct_bins/bin.0.fa"))
                    ])

ref_catdb = channel.fromPath("${params.CAT_ref_db}/${params.cat_db_name}", checkIfExists: true)
ref_cat_diamond = channel.fromPath("${params.CAT_ref_db}/${params.cat_diamond_db_name}", checkIfExists: true)
ref_cat_taxonomy = channel.fromPath("${params.CAT_ref_db}/${params.cat_taxonomy_db}", checkIfExists: true)
ref_gunc = channel.fromPath("${params.gunc_ref_db}", checkIfExists: true)

workflow GGP {
    collect_binners.view()
    BIN_REFINEMENT(collect_binners)
    output_for_prok_part = BIN_REFINEMENT.out.bin_ref_bins.combine(metabat_depth_for_coverage, by: 0)

    CLEAN_AND_FILTER_BINS(output_for_prok_part, ref_catdb.first(), ref_cat_diamond.first(), ref_cat_taxonomy.first(), ref_gunc.first())
}
