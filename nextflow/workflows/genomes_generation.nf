/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/
study_name =            channel.value(params.study_name)
assemblies =            channel.fromPath("${params.assemblies}/*", checkIfExists: true)
raw_reads =             channel.fromPath("${params.raw_reads}/*", checkIfExists: true)
rename_file =           channel.fromPath(params.rename_file)
/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/
ref_genome =            file(params.ref_genome)
ref_genome_index =      file("${params.ref_genome}.*")
ref_eukcc =             channel.fromPath("${params.eukcc_ref_db}", checkIfExists: true)
ref_catdb =             channel.fromPath("${params.CAT_ref_db}/${params.cat_db_name}", checkIfExists: true)
ref_cat_diamond =       channel.fromPath("${params.CAT_ref_db}/${params.cat_diamond_db_name}", checkIfExists: true)
ref_cat_taxonomy =      channel.fromPath("${params.CAT_ref_db}/${params.cat_taxonomy_db}", checkIfExists: true)
ref_gunc =              channel.fromPath("${params.gunc_ref_db}", checkIfExists: true)
ref_checkm =            channel.fromPath("${params.checkm_ref_db}", checkIfExists: true)
ref_rfam_rrna_models =  channel.fromPath("${params.rfam_rrna_models}", checkIfExists: true)
ref_gtdbtk =            channel.fromPath("${params.gtdbtk}", checkIfExists: true)
ref_busco =             channel.fromPath("${params.busco_ref_db}", checkIfExists: true)
/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { PROCESS_INPUT                                    } from '../subworkflows/local/subwf_process_input_files'
include { DECONTAMINATION                                  } from '../subworkflows/local/subwf_decontamination'
include { ALIGN                                            } from '../subworkflows/local/subwf_alignment'
include { QC_AND_MERGE as FILTERING_READS_FASTP            } from '../subworkflows/nf-core/taxprofiler/qc_and_merge'
include { BINNING                                          } from '../subworkflows/nf-core/mag/subwf_mag_binning'
include { EUK_SUBWF                                        } from '../subworkflows/per_sample_euk_part'
include { PROK_SUBWF                                       } from '../subworkflows/per_sample_prok_part'
//include { GZIP                                   } from '../modules/local/utils'
include { GUNZIP_FILES_IN_FOLDER as GUNZIP_BINNER1         } from '../modules/mgbinrefinder/utils'
include { GUNZIP_FILES_IN_FOLDER as GUNZIP_BINNER2         } from '../modules/mgbinrefinder/utils'
include { GUNZIP_FILES_IN_FOLDER as GUNZIP_BINNER3         } from '../modules/mgbinrefinder/utils'
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {
    // ---- combine data for reads and contigs pre-processing
    groupAssemblies = { fasta_file ->
            def cluster = fasta_file.toString().tokenize("/")[-1].tokenize(".")[0]
                def meta = [:]
                meta.id = cluster
                meta.library_layout = "paired"
                meta.single_end = false
            return tuple(meta, fasta_file)
        }
    groupReads = { fastq ->
            def cluster = fastq.toString().tokenize("/")[-1].tokenize(".")[0].tokenize('_')[0]
            def meta = [:]
            meta.id = cluster
            meta.library_layout = "paired"
            meta.single_end = false
            return tuple(meta, fastq)
        }

    tuple_assemblies = assemblies.map(groupAssemblies)      // [ meta, assembly_file ]
    tuple_reads = raw_reads.map(groupReads).groupTuple()    // [ meta, [raw_reads] ]
    data_by_run_accession = tuple_assemblies.combine(tuple_reads, by: 0)  // [ meta, assembly_file, [raw_reads] ]
    data_by_run_accession.view()

    // ---- pre-processing
    PROCESS_INPUT(data_by_run_accession, rename_file)      // output: [ meta, assembly_file, [raw_reads] ]
    assembly = PROCESS_INPUT.out.return_tuple.map(it -> [it[0], it[1]])

    // --- trimming reads
    FILTERING_READS_FASTP(PROCESS_INPUT.out.return_tuple.map(it -> [it[0], it[2]]))

    // --- decontamination
    ref_genome_ch = FILTERING_READS_FASTP.out.reads.map { it -> [it[0], ref_genome, ref_genome_index] }
    DECONTAMINATION(FILTERING_READS_FASTP.out.reads, ref_genome_ch)

    // --- align reads to assembly
    ALIGN(assembly.combine(DECONTAMINATION.out.decontaminated_reads, by:0))  // tuple (meta, fasta, [reads])

    // ---- binning
    BINNING(ALIGN.out.output, DECONTAMINATION.out.decontaminated_reads)
    GUNZIP_BINNER1(BINNING.out.concoct_bins)
    GUNZIP_BINNER2(BINNING.out.maxbin_bins)
    GUNZIP_BINNER3(BINNING.out.metabat_bins)

    if ( !params.skip_euk ) {
        // ---- detect euk
        // input: tuple( meta, assembly_file, [raw_reads], concoct_folder, metabat_folder ), dbs...
        euk_input = assembly.combine(DECONTAMINATION.out.decontaminated_reads, by:0).combine(GUNZIP_BINNER1.out.output_folder, by:0).combine(GUNZIP_BINNER2.out.output_folder, by:0)
        EUK_SUBWF(euk_input, ref_eukcc, ref_busco, ref_catdb, ref_cat_taxonomy)
    }
    if ( !params.skip_prok ) {
        // ---- detect prok
        // input: tuple( meta, concoct, metabat, maxbin, depth_file), dbs...
        prok_input = GUNZIP_BINNER1.out.output_list.combine(GUNZIP_BINNER2.out.output_list, by:0).combine(GUNZIP_BINNER3.out.output_list, by:0).combine(BINNING.out.metabat2depths, by:0)
        PROK_SUBWF(prok_input, ref_catdb, ref_cat_diamond, ref_cat_taxonomy, ref_gunc, ref_checkm, ref_gtdbtk, ref_rfam_rrna_models)
    }
    // ---- compress results
    //GZIP(PROK_SUBWF.out.prok_mags, channel.value("dereplicated_genomes_prok"))
}
