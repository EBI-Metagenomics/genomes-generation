#!/usr/bin/env
cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 1
    ramMin: 200

baseCommand: [ sed ]

inputs:
  input_file:
    type: File
    inputBinding:
      position: 2
  command:
    type: string
    inputBinding:
      position: 1
      # prefix: "-i"
  output_name: string


stdout: $(inputs.output_name)

outputs:
  output_file:
    type: stdout

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