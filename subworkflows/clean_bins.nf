include { CHANGE_DOT_TO_UNDERSCORE } from '../modules/prepare_input'

/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow CLEAN_BINS {
    take:
        mode
        name
        contigs
        reads
        ref_genome
        ref_genome_name
    main:


    emit:
        reads = BEDTOOLS_BAMTOFASTQ.out.reads_cleaned
}
