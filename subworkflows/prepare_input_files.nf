/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP } from '../modules/gunzip'
include { CHANGE_DOT_TO_UNDERSCORE } from '../modules/prepare_data'
include { TRIM_GALORE } from '../modules/prepare_input'
include { MAP_HOST_GENOME } from '../modules/prepare_input'
include { BEDTOOLS_BAMTOFASTQ } from '../modules/bedtools'

workflow PREPARE_INPUT {
    take:
        mode
        name
        contigs
        reads
        ref_genome
        ref_genome_name
    main:
    // fix contig names
    GUNZIP(name, contigs)
    CHANGE_DOT_TO_UNDERSCORE(GUNZIP.out.uncompressed)

    // cleaning reads
    reads_list = reads.collect()

    TRIM_GALORE(mode, name, reads_list)

    MAP_HOST_GENOME(mode, name, ref_genome, ref_genome_name, TRIM_GALORE.out.reads_trimmed)

    BEDTOOLS_BAMTOFASTQ(mode, name, MAP_HOST_GENOME.out.bam_sorted)

    emit:
        contigs_fixed = CHANGE_DOT_TO_UNDERSCORE.out.return_contigs
        reads_cleaned = BEDTOOLS_BAMTOFASTQ.out.reads_cleaned
}
