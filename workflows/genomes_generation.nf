/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/
sample_name = channel.value(params.sample_name)
mode = channel.value(params.mode)
contigs = channel.fromPath(params.contigs, checkIfExists: true)
sr = Channel.fromPath('NO_FILE')
pf = Channel.fromPath('NO_FILE')
pr = Channel.fromPath('NO_FILE')
if (params.mode == 'single') {
    sr = channel.fromPath(params.single_end, checkIfExists: true)
}
else if (params.mode == 'paired') {
    pf = channel.fromPath(params.paired_end_forward, checkIfExists: true)
    pr = channel.fromPath(params.paired_end_reverse, checkIfExists: true)
}

/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/
ref_genome = channel.fromPath(params.ref_genome, checkIfExists: true)
ref_genome_name = channel.value(params.ref_genome_name)

/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { CHANGE_DOT_TO_UNDERSCORE } from '../modules/prepare_input'
include { TRIM_GALORE } from '../modules/prepare_input'
include { MAP_HOST_GENOME } from '../modules/prepare_input'

include { BINNING } from '../modules/metawrap'

/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {

    // CHANGE_DOT_TO_UNDERSCORE(contigs)

    //TRIM_GALORE(mode, sample_name, sr, pf, pr)

    MAP_HOST_GENOME(mode, sample_name, ref_genome, ref_genome_name, sr, pf, pr)
        //TRIM_GALORE.out.single_trimmed,
        //TRIM_GALORE.out.paired_forward_trimmed,
        //TRIM_GALORE.out.paired_reverse_trimmed)

    //BINNING(mode,
    //    contigs,
    //    channel.value(""),
    //    fr, rr
    //)

}
