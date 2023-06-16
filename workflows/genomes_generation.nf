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
include { per_sample_GGP } from '../subworkflows/per_sample_GGP'

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

    // prepare input files
    // input: data_by_run_accession: channel ([ run_accession, assembly_file, [raw_reads] ], ...)
    per_sample_GGP(
        data_by_run_accession,
        ref_genome.first(),
        ref_genome_name,
        rename_file.first(),
        ref_eukcc.first(),
        ref_catdb.first(),
        ref_cat_diamond.first(),
        ref_cat_taxonomy.first(),
        ref_gunc.first(),
        ref_checkm.first(),
        ref_gtdbtk.first(),
        ref_rfam_rrna_models.first()
    )
    // aggregate outputs
    // coverage
    // drep MAGs
    // eukcc MAGs
    // busco MAGs
    // QC MAGs
}
