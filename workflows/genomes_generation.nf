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
ref_eukcc = file("${params.eukcc_ref_db}", checkIfExists: true)
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

    PREPARE_INPUT(data_by_run_accession, ref_genome, ref_genome_name, rename_file)

    BINNING(PREPARE_INPUT.out.contigs_fixed, PREPARE_INPUT.out.reads_cleaned)

    //PROK_SUBWF(sample_name,ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc, ref_checkm, ref_gtdbtk, ref_rfam_rrna_models)

    EUK_SUBWF(PREPARE_INPUT.out.reads_cleaned, PREPARE_INPUT.out.contigs_fixed, BINNING.out.concoct_bins, BINNING.out.metabat2_bins, ref_eukcc.first())
}
