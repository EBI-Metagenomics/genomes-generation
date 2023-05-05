#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GGP } from './workflows/genomes_generation_test'

workflow {
    GGP ()
}
