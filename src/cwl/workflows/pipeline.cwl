#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow


requirements:
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  raw_reads_f:
    type: File[]
  raw_reads_b:
    type: File[]
  scaffolds:
    type: File[]
  input_eukcc_db: Directory
  checkM_db: Directory


outputs: []

steps:

# TODO: cmseq -> coverage
# TODO: remove unbinned files?

  process_sample:
    scatter:
    run: per_sample_pipeline.cwl
    in:
      raw_reads1: raw_reads_f
      raw_reads2: raw_reads_b
      spades_scaffolds: scaffolds
      input_eukcc_db: input_eukcc_db
      checkM_db: checkM_db
    out:
      -



$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20