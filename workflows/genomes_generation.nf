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
include { per_sample_GGP } from '../workflows/per_sample_GGP'

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
    // input: channel ([ run_accession, assembly_file, [raw_reads] ], ...)
    per_sample_GGP(data_by_run_accession, ref_genome.first(), ref_genome_name, rename_file.first(), ref_eukcc.first())
}
