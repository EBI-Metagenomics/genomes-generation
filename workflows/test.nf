/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { BINNING } from '../subworkflows/subwf_mag_binning'

/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {

    //assemblies           // channel: [ val(meta), path(assembly), path(bams), path(bais) ]
    //reads                // channel: [ val(meta), [ reads ] ]

    INPUT_CHECK (
        file(params.input)
    )
    BINNING(INPUT_CHECK.out.assemblies, INPUT_CHECK.out.reads)
}
