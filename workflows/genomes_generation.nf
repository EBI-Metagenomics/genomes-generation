/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/
study_name = channel.value(params.study_name)
assemblies = channel.fromPath("${params.assemblies}/*", checkIfExists: true)
raw_reads = channel.fromPath("${params.raw_reads}/*", checkIfExists: true)
rename_file = channel.fromPath(params.rename_file)
/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/
ref_genome = channel.fromPath(params.ref_genome, checkIfExists: true)
ref_genome_name = channel.value(params.ref_genome_name)
ref_eukcc = channel.fromPath(params.eukcc_ref_db, checkIfExists: true)
ref_catdb = channel.fromPath("${params.CAT_ref_db}/${params.cat_db_name}", checkIfExists: true)
ref_cat_diamond = channel.fromPath("${params.CAT_ref_db}/${params.cat_diamond_db_name}", checkIfExists: true)
ref_cat_taxonomy = channel.fromPath("${params.CAT_ref_db}/${params.cat_taxonomy_db}", checkIfExists: true)
ref_gunc = channel.fromPath("${params.gunc_ref_db}", checkIfExists: true)
ref_checkm = channel.fromPath("${params.checkm_ref_db}", checkIfExists: true)
ref_rfam_rrna_models = channel.fromPath("${params.rfam_rrna_models}", checkIfExists: true)
ref_gtdbtk = channel.fromPath("${params.gtdbtk}", checkIfExists: true)
/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_INPUT } from '../subworkflows/subwf_prepare_input_files'
include { BINNING } from '../subworkflows/subwf_binning'
include { EUK_SUBWF } from '../subworkflows/per_sample_euk_part'
include { PROK_SUBWF } from '../subworkflows/per_sample_prok_part'
include { GZIP } from '../modules/utils'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {
    // ---- combine data for reads and contigs pre-processing
    groupAssemblies = { fasta_file ->
            def cluster = fasta_file.toString().tokenize("/")[-1].tokenize(".")[0]
            return tuple(cluster, fasta_file)
        }
    groupReads = { fastq ->
            def cluster = fastq.toString().tokenize("/")[-1].tokenize(".")[0].tokenize('_')[0]
            return tuple(cluster, fastq)
        }

    tuple_assemblies = assemblies.map(groupAssemblies)      // [ run_accession, assembly_file ]
    tuple_reads = raw_reads.map(groupReads).groupTuple()    // [ run_accession, [raw_reads] ]
    data_by_run_accession = tuple_assemblies.combine(tuple_reads, by: 0)  // [ run_accession, assembly_file, [raw_reads] ]
    data_by_run_accession.view()

    // ---- pre-processing
    PREPARE_INPUT(data_by_run_accession, ref_genome, ref_genome_name, rename_file)      // output: [ run_accession, assembly_file, [raw_reads] ]

    // ---- binning
    BINNING(PREPARE_INPUT.out.return_tuple)

    // ---- detect euk
    // input: tuple( run_accession, assembly_file, [raw_reads], concoct_folder, metabat_folder )
    EUK_SUBWF(BINNING.out.output_for_euk_part, ref_eukcc.first())

    // ---- detect prok
    // input: tuple( run_accession, bin_refinement, depth_file )
    PROK_SUBWF(BINNING.out.output_for_prok_part, ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc, ref_checkm, ref_gtdbtk, ref_rfam_rrna_models)

    //GZIP(PROK_SUBWF.out.prok_mags, channel.value("dereplicated_genomes_prok"))
}
