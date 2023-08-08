//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    //checked_csv = SAMPLESHEET_CHECK(samplesheet).out.csv

    Channel
        .fromPath( samplesheet )
        .splitCsv( header: true, sep: ',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    Channel
        .fromPath( samplesheet )
        .splitCsv( header: true, sep: ',' )
        .map { create_assembly_channel(it) }
        .set { assemblies }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    assemblies                                // channel: [ val(meta), assembly, bam, bai ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.study_name
    //meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

// Function to get list of [ meta, assembly, bam ,bai ]
def create_assembly_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.study_name
    //meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def assembly_meta = []
    assembly_meta = [ meta, file(row.assembly), file(row.bam), file(row.bai) ]
    return assembly_meta
}