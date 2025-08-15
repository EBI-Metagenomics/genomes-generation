include { samplesheetToList                                   } from 'plugin/nf-schema'
include { GENERATE_INPUT_SAMPLESHEET                          } from '../../modules/local/generate_input_samplesheet'
include { DOWNLOAD_FROM_FIRE as DOWNLOAD_FROM_FIRE_READS      } from '../../modules/local/download_from_fire'
include { DOWNLOAD_FROM_FIRE as DOWNLOAD_FROM_FIRE_ASSEMBLIES } from '../../modules/local/download_from_fire'


workflow SAMPLESHEET_GENERATION 
{

    main:

    ch_versions = Channel.empty()

    if (params.samplesheet) {
        samplesheet_ch = Channel.fromList(
           samplesheetToList(
              params.samplesheet,
              "${projectDir}/assets/schema_input.json"
           )
        )
    }
    else {
        GENERATE_INPUT_SAMPLESHEET(
            params.ena_assembly_study_accession,
            params.ena_raw_reads_study_accession
        )
        ch_versions = ch_versions.mix( GENERATE_INPUT_SAMPLESHEET.out.versions )

        samplesheet_ch = GENERATE_INPUT_SAMPLESHEET.out.samplesheet
        .flatMap { samplesheet_file ->
            samplesheetToList(
                samplesheet_file,
                "${projectDir}/assets/schema_input.json"
            )
        }
    }

    // ---- download data from s3 fire for assemblies and reads
    if (params.download_data) {
        assembly_and_runs = samplesheet_ch.map{ meta, fq1, fq2, assembly, _concoct, _metabat, _maxbin, _depth ->
            def is_single_end = (fq2 == [])
            def reads = is_single_end ? [fq1] : [fq1, fq2]
            return tuple(meta + [single_end: is_single_end], assembly, reads)
        }
        DOWNLOAD_FROM_FIRE_READS(
            assembly_and_runs.map{ meta, _assemblies, reads -> [meta, reads] }
            )
        ch_versions = ch_versions.mix( DOWNLOAD_FROM_FIRE_READS.out.versions )

        DOWNLOAD_FROM_FIRE_ASSEMBLIES(
            assembly_and_runs.map{ meta, assemblies, _reads -> [meta, [assemblies]] }
        )
        ch_versions = ch_versions.mix( DOWNLOAD_FROM_FIRE_ASSEMBLIES.out.versions )

        assembly_and_reads = DOWNLOAD_FROM_FIRE_ASSEMBLIES.out.downloaded_files.join( DOWNLOAD_FROM_FIRE_READS.out.downloaded_files )
    } else {
        assembly_and_reads = samplesheet_ch.map{ meta, fq1, fq2, assembly, _concoct, _metabat, _maxbin, _depth ->
            def is_single_end = (fq2 == [])
            def reads = is_single_end ? [fq1] : [fq1, fq2]
            return tuple(meta + [single_end: is_single_end], assembly, reads)
        }
    }

    // --- process samplesheet
    samplesheet_ch.multiMap{ meta, _fq1, fq2, _assembly, concoct, metabat, maxbin, depth ->
        def is_single_end = (fq2 == [])
        concoct: concoct ? tuple(meta + [single_end: is_single_end], concoct) : null
        metabat: metabat ? tuple(meta + [single_end: is_single_end], metabat) : null
        maxbin: maxbin ? tuple(meta + [single_end: is_single_end], maxbin) : null
        jgi_depth: depth ? tuple(meta + [single_end: is_single_end], depth) : null
    }.set {
        input
    }


    // --- assembly software file ---
    if (params.assembly_software_file) {
        assembly_software = file(params.assembly_software_file, checkIfExists: true)
    }
    else {
        assembly_software = assembly_and_reads
        .collectFile(name: 'assembly_software.tsv', newLine: true) { row ->
            def meta = row[0]
            "${meta.id}\t${meta.assembler}"
        }
    }


    emit:
    assembly_and_reads   = assembly_and_reads // channel: [ val(meta), assembly, [ reads ] ]
    concoct_bins         = input.concoct
    metabat_bins         = input.metabat
    maxbin_bins          = input.maxbin
    jgi_depth            = input.jgi_depth
    assembly_software    = assembly_software  // file(id \t assembler)
    versions             = ch_versions        // channel: [ versions.yml ]
}