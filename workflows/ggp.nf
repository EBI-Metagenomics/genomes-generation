/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap    } from 'plugin/nf-schema'  // for multiqc

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { FINALIZE_LOGGING            } from '../modules/local/utils'
include { GUNZIP as GUNZIP_ASSEMBLY   } from '../modules/local/utils'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BINNING              } from '../subworkflows/local/binning'
include { DECONTAMINATION      } from '../subworkflows/local/decontamination'
//include { EUKCC_MERGE          } from '../subworkflows/local/eukcc_merge'
//include { EUK_MAGS_GENERATION  } from '../subworkflows/local/euk_mags_generation'
include { INPUT_PREPROCESSING  } from '../subworkflows/local/input_preprocessing'
//include { PREPARE_UPLOAD_FILES } from '../subworkflows/local/prepare_upload'
//include { PROK_MAGS_GENERATION } from '../subworkflows/local/prok_mags_generation'
include { QC_AND_MERGE_READS   } from '../subworkflows/local/qc_and_merge'

/*
~~~~~~~~~~~~~~~~~~
    DBs
~~~~~~~~~~~~~~~~~~
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

assembly_software  = file(params.assembly_software_file, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Run workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GGP {

    take:
    samplesheet

    main:

    ch_versions    = Channel.empty()
    ch_log         = Channel.empty()

    euk_genomes    = Channel.empty()
    stats_euks     = Channel.empty()
    coverage_euks  = Channel.empty()
    prok_genomes   = Channel.empty()
    stats_proks    = Channel.empty()
    coverage_proks = Channel.empty()
    rna            = Channel.empty()
    taxonomy_euks  = Channel.empty()
    taxonomy_proks = Channel.empty()

    /*
    * ---- combine data for reads, contigs and bins ----
    */
    samplesheet.multiMap{ meta, assembly, fq1, fq2, concoctf, concoctd, metabatf, metabatd, maxbinf, maxbind ->
        def is_single_end = (fq2 == [])
        def reads = is_single_end ? [fq1] : [fq1, fq2]
        assembly_and_runs : tuple(meta + [single_end: is_single_end], assembly, reads)
        concoct: concoctf && concoctd ? tuple(meta + [single_end: is_single_end], concoctf, concoctd, 'concoct') : null
        metabat: metabatf && metabatd ? tuple(meta + [single_end: is_single_end], metabatf, metabatd, 'metabat') : null
        maxbin: maxbinf && maxbind ? tuple(meta + [single_end: is_single_end], maxbinf, maxbind, 'maxbin') : null
    }.set {
        input
    }

    /*
    * --- pre-processing input files ---
    * skip that step with --skip_preprocessing_input
    * change ERR to ERZ in reads
    * change . to _ in assembly files
    */
    INPUT_PREPROCESSING(
        input.assembly_and_runs
    )
    ch_versions = ch_versions.mix(INPUT_PREPROCESSING.out.versions)

    /*
    * --- trimming reads ----
    * merge step is regulated with --merge_pairs (default: false)
    */
    QC_AND_MERGE_READS(
        INPUT_PREPROCESSING.out.assembly_and_reads.map{ meta, _, reads -> [meta, reads] }
    )
    ch_versions = ch_versions.mix( QC_AND_MERGE_READS.out.versions )

    /*
    * --- reads decontamination ----
    * skip that step with --skip_decontamination
    */
    reference_genome         = file(params.ref_genome, checkIfExists: true)
    reference_genome_index   = file("${reference_genome.parent}/*.fa*.*")
    DECONTAMINATION(
        QC_AND_MERGE_READS.out.reads,
        reference_genome,
        reference_genome_index
    )
    ch_versions = ch_versions.mix( DECONTAMINATION.out.versions )

    // --- check filtered reads quality --- //
    FASTQC (
        DECONTAMINATION.out.decontaminated_reads
    )

    /*
    * --- binning ---
    */
    BINNING(
        INPUT_PREPROCESSING.out.assembly_and_reads.map{ meta, assembly, _ -> [meta, assembly] }.join(DECONTAMINATION.out.decontaminated_reads)
    )

    // for multiqc
    // TODO return fastqc from INPUT_PREPROCESSING
    // TODO return fastqc result after decontamination
    // TODO return ALIGN samtools stats

}
