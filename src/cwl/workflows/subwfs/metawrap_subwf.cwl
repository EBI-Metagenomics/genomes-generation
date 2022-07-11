#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

doc: |
      gunzip + metaWRAP


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
  fasta:
    type: File
    outputSource: change_dots_to_underscores/output_file
  bins_concoct:
    type: Directory?
    outputSource: metawrap/concoct_bins
  bins_metabat2:
    type: Directory?
    outputSource: metawrap/metabat2_bins

steps:

  uncompress_fq1:
    run: ../../utils/gunzip.cwl
    in:
      compressed_file: raw_reads1
    out:
      - uncompressed_file

  uncompress_fq2:
    run: ../../utils/gunzip.cwl
    in:
      compressed_file: raw_reads2
    out:
      - uncompressed_file

  change_dots_to_underscores:
    run: ../../utils/sed.cwl
    in:
      input_file: spades_scaffolds
      command: { default: 's/\./_/' }
      output_name:
        source: spades_scaffolds
        valueFrom: $(self.nameroot)_renamed.fasta
    out:
      - output_file

  metawrap:
    run: ../../tools/metawrap/metawrap.cwl
    in:
      contigs: change_dots_to_underscores/output_file
      reads:
        - uncompress_fq1/uncompressed_file
        - uncompress_fq2/uncompressed_file
    out:
      - concoct_bins
      - metabat2_bins


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20