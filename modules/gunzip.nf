/*
    ~~~~~~~~~~~~~~~~~~
     uncompress
    ~~~~~~~~~~~~~~~~~~
*/
process GUNZIP {
    input:
    val name
    path compressed_file

    output:
    path "out/*", emit: uncompressed

    script:
    """
    mkdir out
    cp ${compressed_file} out
    cd out
    gunzip *
    """
}