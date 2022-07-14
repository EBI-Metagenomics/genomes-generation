#!/usr/bin/env
cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 1
    ramMin: 200
  InitialWorkDirRequirement:
    listing:
      - class: File
        location: choose_genomes.sh

inputs:
  concoct_csv:
    type: File
    inputBinding:
      position: 1
      prefix: '-c'

  metabat2_csv:
    type: File
    inputBinding:
      position: 2
      prefix: '-m'

  concoct_merged:
    type: Directory
    inputBinding:
      position: 3
      prefix: '-a'

  metabat2_merged:
    type: Directory
    inputBinding:
      position: 4
      prefix: '-b'

baseCommand: [ bash, choose_genomes.sh ]

outputs:
  input_genomes:
    type: File[]  #Directory
    outputBinding:
      glob: "genomes/*"
  genomes_list:
    type: File
    format: edam:format_1964
    outputBinding:
      glob: "genomes.txt"
  quality:
    type: File
    format: format_3752
    outputBinding:
      glob: "quality.csv"

#hints:
#  - class: DockerRequirement
#    dockerPull: alpine:3.7


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-29