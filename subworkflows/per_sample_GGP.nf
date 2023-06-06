include { PREPARE_INPUT } from '../subworkflows/subwf_prepare_input_files'
include { BINNING } from '../subworkflows/subwf_binning'
include { EUK_SUBWF } from '../subworkflows/per_sample_euk_part'
include { PROK_SUBWF } from '../subworkflows/per_sample_prok_part'

/*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~
     Run per-sample workflow
    ~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow per_sample_GGP {
    take:
        data_by_run_accession   // input: one item of [ run_accession, assembly_file, [raw_reads] ]
        ref_genome
        ref_genome_name
        rename_file
        ref_eukcc
        ref_catdb
        ref_cat_diamond
        ref_cat_taxonomy
        ref_gunc
        ref_checkm
        ref_gtdbtk
        ref_rfam_rrna_models
    main:
    // ---- pre-processing
    PREPARE_INPUT(data_by_run_accession, ref_genome, ref_genome_name, rename_file)      // output: [ run_accession, assembly_file, [raw_reads] ]

    // ---- binning
    BINNING(PREPARE_INPUT.out.return_tuple)

    // ---- detect euk
    // input: tuple( run_accession, assembly_file, [raw_reads], concoct_folder, metabat_folder )
    // done before DREP
    EUK_SUBWF(BINNING.out.output_for_euk_part, ref_eukcc.first())

    // ---- detect prok
    // input: tuple( run_accession, bin_refinement, depth_file )
    PROK_SUBWF(BINNING.out.output_for_prok_part, ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc, ref_checkm, ref_gtdbtk, ref_rfam_rrna_models)
}
