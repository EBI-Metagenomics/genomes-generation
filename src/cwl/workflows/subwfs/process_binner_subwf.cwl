#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}


inputs:
  input_bam:
    type: File
    secondaryFiles: .bai
  input_bins: Directory
  binner_outdir: string
  eukcc_db: Directory

outputs:
  linktable_file:
    type: File?
    outputSource: linktable_step/link_table

  eukcc_out:
    type: Directory
    outputSource: eukCC/eukcc_dir
  eukcc_csv:
    type: File
    outputSource: eukCC/csv

steps:

  linktable_step:
    run: ../../tools/linktable/linktable.cwl
    in:
      bam: input_bam
      bins: input_bins
    out:
      - link_table

  eukCC:
    run: ../../tools/eukcc2/eukcc_binmerging.cwl
    in:
      links: linktable_step/link_table
      outdir: binner_outdir
      db: eukcc_db
      bins: input_bins
    out:
      - eukcc_dir
      - csv


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20