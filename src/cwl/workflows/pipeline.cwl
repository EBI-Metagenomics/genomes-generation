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
  raw_reads1:
    type: File
    format: edam:format_1930  # FASTQ
  raw_reads2:
    type: File
    format: edam:format_1930  # FASTQ
  spades_scaffolds:
    type: File
    format: edam:format_1929 # FASTA
  input_eukcc_db: Directory


outputs:

  bam_out:
    type: File
    outputSource: minimap_align/bamfile

  concoct_binning:
    type: Directory?
    outputSource: binning/bins_concoct
  metabat2_binning:
    type: Directory?
    outputSource: binning/bins_metabat2

  concoct_eukcc:
    type: Directory
    outputSource: process_concoct/eukcc_out
  metabat2_eukcc:
    type: Directory
    outputSource: process_metabat2/eukcc_out

  concoct_linktable:
    type: File?
    outputSource: process_concoct/linktable_file
  metabat2_linktable:
    type: File?
    outputSource: process_metabat2/linktable_file


steps:

  unzip_scaffolds:
    run: ../utils/gunzip.cwl
    in:
      compressed_file: spades_scaffolds
    out:
      - uncompressed_file

  minimap_align:
    run: ../tools/minimap2/minimap2.cwl
    in:
      reads1: raw_reads1
      reads2: raw_reads2
      scaffolds: unzip_scaffolds/uncompressed_file
    out:
      - bamfile

  binning:
    run: subwfs/metawrap_subwf.cwl
    in:
      raw_reads1: raw_reads1
      raw_reads2: raw_reads2
      spades_scaffolds: unzip_scaffolds/uncompressed_file
    out:
      - bins_concoct
      - bins_metabat2

  process_concoct:
    run: subwfs/process_binner_subwf.cwl
    in:
      input_bam: minimap_align/bamfile
      input_bins: binning/bins_concoct
      eukcc_db: input_eukcc_db
      binner_outdir: {default: "concoct_eukcc"}
    out:
      - linktable_file
      - eukcc_out

  process_metabat2:
    run: subwfs/process_binner_subwf.cwl
    in:
      input_bam: minimap_align/bamfile
      input_bins: binning/bins_metabat2
      eukcc_db: input_eukcc_db
      binner_outdir: {default: "metabat2_eukcc"}
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