#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GGP } from './nextflow/workflows/genomes_generation'

workflow {
    GGP ()
}
