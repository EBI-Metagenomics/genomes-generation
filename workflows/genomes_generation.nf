include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}

// Validate input parameters
// validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

if (params.help) {
   log.info paramsHelp("nextflow run ebi-metagenomics/genomes-generation --help")
   exit 0
}

// /*
//     ~~~~~~~
//      Input
//     ~~~~~~~
// */

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
    groupReads = { meta, assembly, fq1, fq2 ->
        if (fq2 == []) {
            return tuple(meta, assembly, [fq1])
        }
        else {
            return tuple(meta, assembly, [fq1, fq2])
        }
    }
    assembly_and_runs = Channel.fromSamplesheet("samplesheet", header: true, sep: ',').map(groupReads) // [ meta, assembly_file, [raw_reads] ]
    assembly_and_runs.view()
    tuple_assemblies = assembly_and_runs.map{ meta, assembly, _ -> tuple(meta, assembly)}

    // ---- pre-processing ---- //
    PROCESS_INPUT( assembly_and_runs ) // output: [ meta, assembly_file, [raw_reads] ]

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

        EUK_MAGS_GENERATION( 
            euk_input,
            eukcc_db,
            busco_db,
            cat_db_folder,
            cat_diamond_db,
            cat_taxonomy_db
        )
    }

    if ( !params.skip_prok ) {

        // input: tuple( meta, concoct, metabat, maxbin, depth_file)
        collected_binners_and_depth = concoct_collected_bins.join( maxbin_collected_bins ) \
            .join( metabat_collected_bins ) \
            .join( BINNING.out.metabat2depths )

        PROK_MAGS_GENERATION(
            collected_binners_and_depth,
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
