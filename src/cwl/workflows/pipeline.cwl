#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

doc: |


requirements:
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  raw_reads1:
    type: File
    format: edam:format_1930  # FASTQ
  raw_reads2:
    type: File
    format: edam:format_1930  # FASTQ
  spades_scaffolds:
    type: File
    format: edam:format_1929 # FASTA


outputs:

  checkm_stderr:
    type: File?
    outputSource: checkm_subwf/checkm_err

steps:

  minimap_align:
    run: ../tools/minimap/minimap.cwl
    in:
      reads1: raw_reads1
      reads2: raw_reads2
      scaffolds: spades_scaffolds
    out:
      - bamfile

  binning:
    run: subwfs/metawrap_subwf.cwl
    in:
      raw_reads1: raw_reads1
      raw_reads2: raw_reads2
      spades_scaffolds: spades_scaffolds
    out:
      - bins_concoct
      - bins_metabat2

  process_concoct:
    run: subwfs/process_binner_subwf.cwl
    in:

    out:

  process_metabat2:
    run: subwfs/process_binner_subwf.cwl
    in:

    out:


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20