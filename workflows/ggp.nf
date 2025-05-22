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
//include { MULTIQC                     } from '../modules/nf-core/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BINNING              } from '../subworkflows/local/binning'
include { DECONTAMINATION      } from '../subworkflows/local/decontamination'
//include { EUKCC_MERGE          } from '../subworkflows/local/eukcc_merge'
include { EUK_MAGS_GENERATION  } from '../subworkflows/local/eukaryotic_mags_generation'
include { INPUT_PREPROCESSING  } from '../subworkflows/local/input_preprocessing'
//include { PREPARE_UPLOAD_FILES } from '../subworkflows/local/prepare_upload'
//include { PROK_MAGS_GENERATION } from '../subworkflows/local/prok_mags_generation'
include { QC_AND_MERGE_READS   } from '../subworkflows/local/qc_and_merge'

/*
~~~~~~~~~~~~~~~~~~
    DBs
~~~~~~~~~~~~~~~~~~
*/
//cat_db_folder      = file(params.cat_db_folder, checkIfExists: true)
//cat_diamond_db     = file(params.cat_diamond_db, checkIfExists: true)
//cat_taxonomy_db    = file(params.cat_taxonomy_db, checkIfExists: true)
//gunc_db            = file(params.gunc_db, checkIfExists: true)
//checkm2_db         = file(params.checkm2_db, checkIfExists: true)
//busco_db           = file(params.busco_db, checkIfExists: true)
//gtdbtk_db          = file(params.gtdbtk_db, checkIfExists: true)
//rfam_rrna_models   = file(params.rfam_rrna_models, checkIfExists: true)

//assembly_software  = file(params.assembly_software_file, checkIfExists: true)

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
    //ch_log         = Channel.empty()

    //euk_genomes    = Channel.empty()
    //stats_euks     = Channel.empty()
    //coverage_euks  = Channel.empty()
    //prok_genomes   = Channel.empty()
    //stats_proks    = Channel.empty()
    //coverage_proks = Channel.empty()
    //rna            = Channel.empty()
    //taxonomy_euks  = Channel.empty()
    //taxonomy_proks = Channel.empty()

    /*
    * ---- combine data for reads, contigs and bins ----
    */
    samplesheet.multiMap{ meta, assembly, fq1, fq2, concoctf, concoctd, metabatf, metabatd, maxbinf, maxbind ->
        def is_single_end = (fq2 == [])
        def reads = is_single_end ? [fq1] : [fq1, fq2]
        assembly_and_reads : tuple(meta + [single_end: is_single_end], assembly, reads)
        concoct: concoctf && concoctd ? tuple(meta + [single_end: is_single_end], concoctf, concoctd) : null
        metabat: metabatf && metabatd ? tuple(meta + [single_end: is_single_end], metabatf, metabatd) : null
        maxbin: maxbinf && maxbind ? tuple(meta + [single_end: is_single_end], maxbinf, maxbind) : null
    }.set {
        input
    }
    concoct_sample_ids = input.concoct.filter { it != null }
        .map { meta, f, d -> meta.id }
        .collect()
        .ifEmpty([])
    metabat_sample_ids = input.metabat.filter { it != null }
        .map { meta, f, d -> meta.id }
        .collect()
        .ifEmpty([])
    maxbin_sample_ids = input.maxbin.filter { it != null }
        .map { meta, f, d -> meta.id }
        .collect()
        .ifEmpty([])

    tool_availability = input.assembly_and_reads
        .combine(concoct_sample_ids.map { [it] })
        .combine(metabat_sample_ids.map { [it] })
        .combine(maxbin_sample_ids.map { [it] })
        .map { meta, assembly, reads, concoct_list, metabat_list, maxbin_list ->
            def tools_to_run = []
            if (concoct_list.contains(meta.id)) tools_to_run.add('concoct')
            if (metabat_list.contains(meta.id)) tools_to_run.add('metabat')
            if (maxbin_list.contains(meta.id)) tools_to_run.add('maxbin')
            return [meta, assembly, reads, tools_to_run]
        }

    tool_availability
        .multiMap { meta, assembly, reads, tools ->
            run_concoct: !tools.contains('concoct') ? [meta, assembly, reads] : null
            run_metabat: !tools.contains('metabat') ? [meta, assembly, reads] : null  
            run_maxbin: !tools.contains('maxbin') ? [meta, assembly, reads] : null
        }
        .set { branched_for_tools }

    //branched_for_tools.run_maxbin.filter { it != null }

    BINNING(
        branched_for_tools.run_concoct,
        branched_for_tools.run_metabat,
        branched_for_tools.run_maxbin,
    )

}
