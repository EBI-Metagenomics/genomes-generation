include { samplesheetToList                   } from 'plugin/nf-schema'
include { GENERATE_INPUT_SAMPLESHEET          } from '../../modules/local/generate_input_samplesheet'


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
    input_data        = assembly_and_runs  // channel: [ val(meta), assembly, [ reads ] ]
    assembly_software = assembly_software  // file(id \t assembler)
    versions          = ch_versions        // channel: [ versions.yml ]
}