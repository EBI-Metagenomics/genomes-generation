#!/usr/bin/env
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 1
    ramMin: 200

inputs:
  compressed_file:
    type: File
    inputBinding:
      position: 1
      prefix: '-c'

baseCommand: [ gunzip ]

stdout: $(inputs.compressed_file.nameroot)

outputs:
  uncompressed_file:
    type: stdout
    format: edam:format_1930 # FASTQ

hints:
  - class: DockerRequirement
    dockerPull: alpine:3.7


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-29