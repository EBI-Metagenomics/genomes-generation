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
  raw_reads2:
    type: File
  spades_scaffolds:
    type: File
  input_eukcc_db: Directory
  checkM_db: Directory


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

  binrefine_bins_out:
    type: Directory
    outputSource: binning/binrefine_bins

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

# TODO: cmseq -> coverage
# TODO: remove unbinned files?

  unzip_scaffolds:
    run: ../utils/gunzip.cwl
    in:
      compressed_file: spades_scaffolds
    out:
      - uncompressed_file

  change_dots_to_underscores:
    run: ../utils/sed.cwl
    in:
      input_file: unzip_scaffolds/uncompressed_file
      command: { default: 's/\./_/' }
      output_name:
        source: spades_scaffolds
        valueFrom: $(self.nameroot)_renamed.fasta
    out:
      - output_file

  minimap_align:
    run: ../tools/minimap2/minimap2.cwl
    in:
      reads1: raw_reads1
      reads2: raw_reads2
      scaffolds: change_dots_to_underscores/output_file
    out:
      - bamfile

  binning:
    run: subwfs/metawrap_subwf.cwl
    in:
      raw_reads1: raw_reads1
      raw_reads2: raw_reads2
      spades_scaffolds: change_dots_to_underscores/output_file
      checkM_db: checkM_db
    out:
      - bins_concoct
      - bins_metabat2
      - binrefine_bins

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
      - eukcc_csv

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
      - eukcc_csv

$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20