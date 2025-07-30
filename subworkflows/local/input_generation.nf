include { samplesheetToList                                   } from 'plugin/nf-schema'
include { GENERATE_INPUT_SAMPLESHEET                          } from '../../modules/local/generate_input_samplesheet'
include { DOWNLOAD_FROM_FIRE as DOWNLOAD_FROM_FIRE_READS      } from '../../modules/local/download_from_fire'
include { DOWNLOAD_FROM_FIRE as DOWNLOAD_FROM_FIRE_ASSEMBLIES } from '../../modules/local/download_from_fire'


workflow INPUT_GENERATION {

    ch_versions = Channel.empty()

    // --- getting assemblies anf runs ---
    groupReads = { meta, fq1, fq2, assembly ->
        if (fq2 == []) {
            return tuple(meta + [single_end: true], assembly, [fq1])
        }
        else {
            return tuple(meta + [single_end: false], assembly, [fq1, fq2])
        }
    }

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

    assembly_and_runs = samplesheet_ch.map(groupReads)

    // ---- download data from s3 fire
    if (params.download_data) {
        DOWNLOAD_FROM_FIRE_READS(assembly_and_runs.map{meta, assemblies, reads -> [meta, reads]})
        ch_versions = ch_versions.mix( DOWNLOAD_FROM_FIRE_READS.out.versions )
        DOWNLOAD_FROM_FIRE_ASSEMBLIES(assembly_and_runs.map{meta, assemblies, reads -> [meta, [assemblies]]})
        ch_versions = ch_versions.mix( DOWNLOAD_FROM_FIRE_ASSEMBLIES.out.versions )
        input_data = DOWNLOAD_FROM_FIRE_ASSEMBLIES.out.downloaded_files.join(DOWNLOAD_FROM_FIRE_READS.out.downloaded_files)
    } else {
        input_data = assembly_and_runs
    }


    // --- assembly software file ---
    if (params.assembly_software_file) {
        assembly_software = file(params.assembly_software_file, checkIfExists: true)
    }
    else {
        assembly_software = assembly_and_runs
        .collectFile(name: 'assembly_software.tsv', newLine: true) { row ->
            def meta = row[0]
            "${meta.id}\t${meta.assembler}"
        }
    }

    emit:
    input_data        = input_data         // channel: [ val(meta), assembly, [ reads ] ]
    assembly_software = assembly_software  // file(id \t assembler)
    versions          = ch_versions        // channel: [ versions.yml ]
}