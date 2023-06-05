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
/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_INPUT } from '../subworkflows/prepare_input_files'
include { BINNING } from '../subworkflows/binning'
include { EUK_SUBWF } from './euk_part'
//include { PROK_SUBWF } from './prok_part'

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
    PREPARE_INPUT(data_by_run_accession, ref_genome, ref_genome_name, rename_file)

    // ---- binning
    prepared_data_by_run_accession = PREPARE_INPUT.out.reads_cleaned.combine(PREPARE_INPUT.out.contigs_fixed, by: 0)
    prepared_data_by_run_accession.view()
    // tuple(name, [reads], contigs)
    BINNING(prepared_data_by_run_accession)

    //PROK_SUBWF(sample_name,ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc, ref_checkm, ref_gtdbtk, ref_rfam_rrna_models)

    // ---- combine data for EUK PART
    data_after_binning = prepared_data_by_run_accession.combine(BINNING.out.concoct_bins, by:0).combine(BINNING.out.metabat2_bins, by:0)
    reads = data_after_binning.map(item -> tuple(item[0], item[1]))
    contigs = data_after_binning.map(item -> tuple(item[0], item[2]))
    concoct_bins = data_after_binning.map(item -> tuple(item[0], item[3]))
    metabat2_bins = data_after_binning.map(item -> tuple(item[0], item[4]))

    // ---- detect euk
    EUK_SUBWF(reads, contigs, concoct_bins, metabat2_bins, ref_eukcc.first())
}
