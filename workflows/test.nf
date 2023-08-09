/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK_ALIGN } from '../subworkflows/local/input_check'
include { ALIGN_META } from '../subworkflows/subwf_alignment'
include { BINNING } from '../subworkflows/subwf_mag_binning'

/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {

    INPUT_CHECK_ALIGN (
        file(params.input)
    )
    //assemblies           // channel: [ val(meta), path(assembly)]
    //reads                // channel: [ val(meta), [ reads ] ]
    INPUT_CHECK_ALIGN.out.data.view()
    ALIGN_META(INPUT_CHECK_ALIGN.out.data)

    reads = INPUT_CHECK_ALIGN.out.data.map {it -> [it[0], it[2]]}
    BINNING(ALIGN_META.out.bams, reads)
}
