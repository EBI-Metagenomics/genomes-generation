#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GGP } from './workflows/test2'

workflow {
    GGP ()
}
