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
  eukcc_db: Directory


outputs:

  bam_out:
    type: File
    outputSource: minimap_align/bamfile

  concoct_binning:
    type: Directory
    outputSource: binning/bins_concoct
  metabat2_binning:
    type: Directory
    outputSource: binning/bins_metabat2

  concoct_eukcc:
    type: Directory
    outputSource: process_concoct/eukcc_out
  metabat2_eukcc:
    type: Directory
    outputSource: process_metabat2/eukcc_out

  concoct_linktable:
    type: Directory
    outputSource: process_concoct/linktable_file
  metabat2_linktable:
    type: Directory
    outputSource: process_metabat2/linktable_file

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
      bam:
      bins:
      eukcc_db: eukcc_db
    out:
      - linktable_file
      - eukcc_out

  process_metabat2:
    run: subwfs/process_binner_subwf.cwl
    in:
      bam:
      bins:
      eukcc_db: eukcc_db
    out:
      - linktable_file
      - eukcc_out

$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20