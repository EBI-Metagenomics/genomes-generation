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
  raw_reads1: File
  raw_reads2: File
  spades_scaffolds: File
  checkM_db: Directory


outputs:

  bins_concoct:
    type: Directory?
    outputSource: create_dir_concoct/out
  bins_metabat2:
    type: Directory?
    outputSource: create_dir_metabat2/out

  binrefine_bins:
    type: Directory
    outputSource: binrefine/bins

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

# return without unbinned.fa
  metawrap:
    run: ../../tools/metawrap/metawrap.cwl
    in:
      contigs: spades_scaffolds
      reads:
        - uncompress_fq1/uncompressed_file
        - uncompress_fq2/uncompressed_file
    out:
      - concoct_bins
      - metabat2_bins

  create_dir_concoct:
    run: ../../utils/return_directory.cwl
    in:
      list: metawrap/concoct_bins
      dir_name: {default: "concoct_bins"}
    out: [ out ]

  create_dir_metabat2:
    run: ../../utils/return_directory.cwl
    in:
      list: metawrap/metabat2_bins
      dir_name: {default: "metabat2_bins"}
    out: [ out ]

  binrefine:
    run: ../../tools/metawrap/binrefine.cwl
    in:
      concoct_bin_dir: create_dir_concoct/out
      metabat_bin_dir: create_dir_metabat2/out
      refdata: checkM_db
    out:
      - bins


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20