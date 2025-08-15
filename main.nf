#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ebi-metagenomics/genomes-generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ebi-metagenomics/genomes-generation
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FINALIZE_LOGGING            } from './modules/local/utils'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'

include { SAMPLESHEET_GENERATION      } from './subworkflows/local/samplesheet_generation'
include { PIPELINE_INITIALISATION     } from './subworkflows/local/pipeline_initialisation'

include { GGP                         } from './workflows/ggp'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    ch_versions = Channel.empty()

    PIPELINE_INITIALISATION (
        args
    )
    ch_versions = ch_versions.mix( PIPELINE_INITIALISATION.out.versions )

    SAMPLESHEET_GENERATION ()
    ch_versions = ch_versions.mix( SAMPLESHEET_GENERATION.out.versions )

    GGP (
        SAMPLESHEET_GENERATION.out.assembly_and_reads,
        SAMPLESHEET_GENERATION.out.concoct_bins,
        SAMPLESHEET_GENERATION.out.metabat_bins,
        SAMPLESHEET_GENERATION.out.maxbin_bins,
        SAMPLESHEET_GENERATION.out.jgi_depth,
        SAMPLESHEET_GENERATION.out.assembly_software
    )
    ch_versions = ch_versions.mix( GGP.out.versions )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    FINALIZE_LOGGING (
        GGP.out.pipeline_logging, 
        "structured_pipeline_logging_by_runs.txt"
    )
}
