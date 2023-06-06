include { PREPARE_INPUT } from '../subworkflows/prepare_input_files'
include { BINNING } from '../subworkflows/binning'
//include { EUK_SUBWF } from './euk_part'
//include { PROK_SUBWF } from './prok_part'

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
    main:
    // ---- pre-processing
    accession_reads = data_by_run_accession.map(item -> tuple(item[0], item[2]))
    accession_contigs = data_by_run_accession.map(item -> tuple(item[0], item[1]))
    PREPARE_INPUT(accession_reads, accession_contigs, ref_genome, ref_genome_name, rename_file)

    // ---- binning
    BINNING(PREPARE_INPUT.out.reads_cleaned, PREPARE_INPUT.out.contigs_fixed)

    //PROK_SUBWF(sample_name,ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc, ref_checkm, ref_gtdbtk, ref_rfam_rrna_models)

    // ---- combine data for EUK PART
    //data_after_binning = prepared_data_by_run_accession.combine(BINNING.out.concoct_bins, by:0).combine(BINNING.out.metabat2_bins, by:0)
    //reads = data_after_binning.map(item -> tuple(item[0], item[1]))
    //contigs = data_after_binning.map(item -> tuple(item[0], item[2]))
    //concoct_bins = data_after_binning.map(item -> tuple(item[0], item[3]))
    //metabat2_bins = data_after_binning.map(item -> tuple(item[0], item[4]))

    // ---- detect euk
    //EUK_SUBWF(reads, contigs, concoct_bins, metabat2_bins, ref_eukcc)
}
