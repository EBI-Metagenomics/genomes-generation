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

assemblies  = channel.fromPath("${params.assemblies}/*", checkIfExists: true)
raw_reads   = channel.fromPath("${params.raw_reads}/*", checkIfExists: true)
erz_to_err_mapping_file = file(params.erz_to_err_mapping_file, checkIfExists: true)

/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~process {

    cpus   = { check_max( 1    * task.attempt, 'cpus'   ) }
    memory = { check_max( 6.GB * task.attempt, 'memory' ) }
    time   = { check_max( 4.h  * task.attempt, 'time'   ) }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'
*/
eukcc_db           = file(params.eukcc_db, checkIfExists: true)
cat_db_folder      = file(params.cat_db_folder, checkIfExists: true)
cat_diamond_db     = file(params.cat_diamond_db, checkIfExists: true)
cat_taxonomy_db    = file(params.cat_taxonomy_db, checkIfExists: true)
gunc_db            = file(params.gunc_db, checkIfExists: true)
checkm2_db         = file(params.checkm2_db, checkIfExists: true)
busco_db           = file(params.busco_db, checkIfExists: true)
gtdbtk_db          = file(params.gtdbtk_db, checkIfExists: true)
rfam_rrna_models   = file(params.rfam_rrna_models, checkIfExists: true)

ref_genome         = file(params.ref_genome)
ref_genome_index   = file("${ref_genome.parent}/*.fa.*")

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

    // TODO: use a samplesheet instead of this grouping //
    groupAssemblies = { fasta_file ->
        id = fasta_file.toString().tokenize("/")[-1].tokenize(".")[0]
        meta = [id:id, single_end: false]
        return tuple(meta, fasta_file)
    }
    groupReads = { fastq ->
        id = fastq.toString().tokenize("/")[-1].tokenize(".")[0].tokenize('_')[0]
        meta = [id:id, single_end: false]
        return tuple(meta, fastq)
    }

    tuple_assemblies = assemblies.map( groupAssemblies ) // [ meta, assembly_file ]
    tuple_reads = raw_reads.map( groupReads ).groupTuple() // [ meta, [raw_reads] ]

    assembly_and_runs = tuple_assemblies.join( tuple_reads )  // [ meta, assembly_file, [raw_reads] ]

    // ---- pre-processing ---- //
    PROCESS_INPUT( assembly_and_runs, erz_to_err_mapping_file ) // output: [ meta, assembly_file, [raw_reads] ]

    // --- trimming reads ---- //
    QC_AND_MERGE_READS( PROCESS_INPUT.out.assembly_and_reads.map { meta, _, reads -> [meta, reads] } )

    // --- decontamination ---- //
    // We need a tuple as the alignment and decontamination module needs the input like that
    DECONTAMINATION( QC_AND_MERGE_READS.out.reads, ref_genome, ref_genome_index )

    // For paired end reads and with merge_paris in true, we need to
    // switch the assembly single_end -> true to be in sync with how
    // fastp works. 
    if ( params.merge_pairs ) {
        tuple_assemblies = tuple_assemblies.map { meta, assembly ->
            [ meta + [single_end: true], assembly ]
        }
    }

    // --- align reads to assemblies ---- //
    assembly_and_reads = tuple_assemblies.join( DECONTAMINATION.out.decontaminated_reads )

    ALIGN( assembly_and_reads ) // tuple (meta, fasta, [reads])

    // ---- binning ---- //
    BINNING( ALIGN.out.assembly_bam )

    collectBinsFolder = { meta, bin_folder ->
        [ meta, bin_folder.listFiles().flatten() ]
    }
    concoct_collected_bins = BINNING.out.concoct_bins.map( collectBinsFolder )
    metabat_collected_bins = BINNING.out.metabat_bins.map( collectBinsFolder )
    maxbin_collected_bins = BINNING.out.maxbin_bins.map( collectBinsFolder )

    if ( !params.skip_euk ) {

        // ---- detect euk ---- //
        // input: tuple( meta, assembly_file, [raw_reads], concoct_folder, metabat_folder ), dbs...
        euk_input = tuple_assemblies.join(
            DECONTAMINATION.out.decontaminated_reads
        ).join(
            concoct_collected_bins
        ).join(
            metabat_collected_bins
        )

        EUK_MAGS_GENERATION( euk_input, eukcc_db, busco_db, cat_diamond_db, cat_taxonomy_db )
    }

    if ( !params.skip_prok ) {

        // input: tuple( meta, concoct, metabat, maxbin, depth_file), dbs...
        prok_input = concoct_collected_bins.join( maxbin_collected_bins ) \
            .join( metabat_collected_bins ) \
            .join( BINNING.out.metabat2depths )

        PROK_MAGS_GENERATION(
            prok_input,
            cat_db_folder,
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
