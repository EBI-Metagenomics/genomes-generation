/*
    ~~~~~~~~~~~~~~~~~~
     Run subworkflow
    ~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP } from '../modules/utils'
include { CHANGE_DOT_TO_UNDERSCORE as CHANGE_DOT_TO_UNDERSCORE_CONTIGS} from '../modules/utils'
include { CHANGE_DOT_TO_UNDERSCORE_READS } from '../modules/utils'
include { FASTP as FILTERING_READS_FASTP} from '../modules/fastp'
include { DECONTAMINATION } from '../subworkflows/decontamination'
include { CHANGE_ERR_TO_ERZ as CHANGE_ERR_TO_ERZ_READS } from '../modules/utils'

workflow PREPARE_INPUT {
    take:
        input_data
        ref_genome
        ref_genome_name
        rename_file
    main:
    // --- fix contig names
    // GUNZIP(name, contigs)
    CHANGE_DOT_TO_UNDERSCORE_CONTIGS(input_data.map{item -> tuple(item[0], item[1])})           // tuple(accession, assembly)

    // change ERR in reads to ERZ
    CHANGE_ERR_TO_ERZ_READS(input_data.map{item -> tuple(item[0], item[2])}, rename_file.first())    // tuple(accession, [reads]])
    CHANGE_DOT_TO_UNDERSCORE_READS(CHANGE_ERR_TO_ERZ_READS.out.return_files)

    // --- trimming reads
    FILTERING_READS_FASTP(CHANGE_DOT_TO_UNDERSCORE_READS.out.underscore_reads)

    DECONTAMINATION(FILTERING_READS_FASTP.out.output_reads, ref_genome.first(), ref_genome_name)
    emit:
        contigs_fixed = CHANGE_DOT_TO_UNDERSCORE_CONTIGS.out.underscore_contigs
        reads_cleaned = DECONTAMINATION.out.decontaminated_reads
}
