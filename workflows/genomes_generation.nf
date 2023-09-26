include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

if (params.help) {
   log.info paramsHelp("nextflow run ebi-metagenomics/genomes-generation --help")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// /*
//     ~~~~~~~
//      Input
//     ~~~~~~~
// */

assemblies  = channel.fromPath("$assemblies/*", checkIfExists: true)
raw_reads   = channel.fromPath("$raw_reads}/*", checkIfExists: true)
erz_to_err_mapping_file = channel.fromPath("erz_to_err_mapping_file")

/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/
ref_genome         = file(params.ref_genome)
ref_genome_index   = file("${params.ref_genome}.*")

eukcc_db           = file("${params.eukcc_ref_db}", checkIfExists: true)

cat_diamond_db     = file("${params.cat_diamond_db}", checkIfExists: true)
cat_taxonomy_db    = file("${params.cat_taxonomy_db}", checkIfExists: true)

gunc_db            = file("${params.gunc_ref_db}", checkIfExists: true)

checkm2_db          = file("${params.checkm2_ref_db}", checkIfExists: true)

rfam_rrna_models   = file("${params.rfam_rrna_models}", checkIfExists: true)

ref_gtdbtk         = file("${params.gtdbtk}", checkIfExists: true)

ref_busco          = file("${params.busco_ref_db}", checkIfExists: true)
gtdbtk_db          = file("${params.gtdbtk_db}", checkIfExists: true)


/*
    ~~~~~~~~~~~~~~~~~~
     Subworkflows
    ~~~~~~~~~~~~~~~~~~
*/
include { PROCESS_INPUT        } from '../subworkflows/local/process_input_files'
include { DECONTAMINATION      } from '../subworkflows/local/decontamination'
include { ALIGN                } from '../subworkflows/local/alignment'
include { EUK_MAGS_GENERATION  } from '../subworkflows/local/euk_mags_generation'
include { PROK_MAGS_GENERATION } from '../subworkflows/local/prok_mags_generation'

include { QC_AND_MERGE_READS   } from '../subworkflows/local/qc_and_merge'
include { BINNING              } from '../subworkflows/local/mag_binning'

/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {
    // ---- combine data for reads and contigs pre-processing ---- //
    groupAssemblies = { fasta_file ->
        def cluster = fasta_file.toString().tokenize("/")[-1].tokenize(".")[0]
        def meta = [:]
        meta.id = cluster
        return tuple(meta, fasta_file)
    }
    groupReads = { fastq ->
        def cluster = fastq.toString().tokenize("/")[-1].tokenize(".")[0].tokenize('_')[0]
        def meta = [:]
        meta.id = cluster
        return tuple(meta, fastq)
    }

    tuple_assemblies = assemblies.map(groupAssemblies) // [ meta, assembly_file ]
    tuple_reads = raw_reads.map(groupReads).groupTuple() // [ meta, [raw_reads] ]
    data_by_run_accession = tuple_assemblies.combine(tuple_reads, by: [0])  // [ meta, assembly_file, [raw_reads] ]
    data_by_run_accession.view()

    // ---- pre-processing ---- //
    PROCESS_INPUT( data_by_run_accession, erz_to_err_mapping_file ) // output: [ meta, assembly_file, [raw_reads] ]

    // --- trimming reads ---- //
    QC_AND_MERGE_READS( PROCESS_INPUT.out.assembly_and_reads.map { it -> [it[0], it[2]] } )

    // --- decontamination ---- //
    // We need a tuple as the alignment and decontamination module needs the input like that
    DECONTAMINATION( QC_AND_MERGE_READS.out.reads.map { it -> tuple( it[0], it[1], ref_genome, ref_genome_index ) } )

    // --- align reads to assembly ---- //
    ALIGN( assembly.combine( DECONTAMINATION.out.decontaminated_reads, by: [0]) ) // tuple (meta, fasta, [reads])

    // ---- binning ---- //
    BINNING( ALIGN.out.output, DECONTAMINATION.out.decontaminated_reads )

    concoct_bins = BINNING.out.concoct_bins  // uncompressed folders bins
    maxbin_bins = BINNING.out.maxbin_bins
    metabat_bins = BINNING.out.metabat_bins

    if ( !params.skip_euk ) {
        // ---- detect euk ---- //
        // input: tuple( meta, assembly_file, [raw_reads], concoct_folder, metabat_folder ), dbs...
        euk_input = assembly.combine(
            DECONTAMINATION.out.decontaminated_reads, by: [0]
        ).combine(
            concoct_bins, by: [0]
        ).combine(
            metabat_bins, by: [0]
        )

        EUK_MAGS_GENERATION( euk_input, eukcc_db, busco_db, cat_db, cat_taxonomy_db )
    }

    if ( !params.skip_prok ) {
        // ---- detect prok ---- //
        concoct_list = concoct_bins.map{ it -> [it[0], it[1].listFiles().flatten()] }
        maxbin_list = maxbin_bins.map{ it -> [it[0], it[1].listFiles().flatten()] }
        metabat_list = metabat_bins.map{ it -> [it[0], it[1].listFiles().flatten()] }

        // input: tuple( meta, concoct, metabat, maxbin, depth_file), dbs...
        prok_input = concoct_list.combine( maxbin_list, by: [0] ) \
            .combine( metabat_list, by: [0] ) \
            .combine( BINNING.out.metabat2depths, by: [0] )

        PROK_MAGS_GENERATION(
            prok_input,
            cat_db,
            cat_diamond_db,
            cat_taxonomy_db,
            gunc_db,
            checkm2_db,
            gtdbtk_db,
            rfam_rrna_models
        )
    }

    // ---- compress results ---- //
    //GZIP(PROK_SUBWF.out.prok_mags, channel.value("dereplicated_genomes_prok"))
}
