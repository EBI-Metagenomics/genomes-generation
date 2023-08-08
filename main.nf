#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GGP } from './workflows/test.nf'
//'./workflows/genomes_generation'

workflow {
    GGP ()
}
