//
// Change reads and assemblies input files if necessary
//

include { RENAME_CONTIGS          } from '../../modules/local/rename_contigs'

include { FASTQC                  } from '../../modules/nf-core/fastqc/main'
inclide { SEQKIT_PAIR             } from '../../modules/nf-core/seqkit/pair'
inclide { SEQKIT_SANA             } from '../../modules/nf-core/seqkit/sana'

include { FASTQ_TRIM_FASTP_FASTQC } from '../nf-core/fastq_trim_fastp_fastqc'


workflow INPUT_QC {

    take:
    input_data  // tuple( meta, assembly_file, [raw_reads] )

    main:

    ch_versions = Channel.empty()

    reads = input_data.map { meta, _assembly, reads -> [meta, reads] }
    contigs = input_data.map { meta, assembly, _reads_item -> [meta, assembly] }

    //
    // --- check input reads quality ---
    //
    FASTQC(
        reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // --- skip reads and contigs modifications --
    //
    if ( params.skip_preprocessing_input ) {
        println('skipping pre-processing')
        result = input_data
    }
    //
    // --- add correct accessions
    //
    else {
        /*
        * --- MODIFY CONTIGS ( change headers to sorter names including assembly_accession )
        *  output: tuple(meta, assembly)
        */
        RENAME_CONTIGS(
            contigs
        )
        ch_versions = ch_versions.mix(RENAME_CONTIGS.out.versions)

        /*
        * --- CHECK READS with SEQKIT
        *  input: tuple(meta, [reads]])
        *  output: tuple(meta, [modified_pe_reads]/modified_se_reads])
        */
        SEQKIT_SANA(
            reads
        )
        ch_versions = ch_versions.mix(SEQKIT_SANA.out.versions)

        SEQKIT_PAIR(
            SEQKIT_SANA.out.reads_sanitised
        )
        ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions)

        /* --- fastqc_raw - fastp - fastqc_trim
        // fastp args: reads, [], false, params.merge_pairs (should be false in config)
        FASTQ_TRIM_FASTP_FASTQC(
            SEQKIT_PAIR.out.reads,
            [],
            false?,
            false,
            params.merge_pairs ??? - is it correct arg? I need to turn off merging in fastp,
            false,
            false
        )
        ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)

        result = RENAME_CONTIGS.out.contigs_renamed.join(FASTQ_TRIM_FASTP_FASTQC.out.reads)
    }

    emit:
    assembly_and_reads     = result // tuple( meta, assembly_file, [raw_reads] )
    fastqc_raw_gzip        = FASTQC.out.zip                               // initial reads qc
    //fastqc_seqkit_gzip     = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip   // after seqkit sana and pair qc
    //fastqc_trim_gzip       = FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip  // after fastp qc
    versions               = ch_versions
}
