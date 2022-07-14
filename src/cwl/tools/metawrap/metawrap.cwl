#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "metaWrap binning tool (binning)"

hints:
  DockerRequirement:
    dockerPull: "quay.io/microbiome-informatics/metawrap:latest"

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: $(inputs.reads)
  ResourceRequirement:
    coresMin: 4

baseCommand: [ 'metawrap', 'binning', '--concoct', '--metabat2']

arguments:
  - valueFrom: $(runtime.ram)
    prefix: -m
  - valueFrom: $(runtime.cores)
    prefix: -t

inputs:
  contigs:
    type: File
    inputBinding:
      prefix: "-a"
  reads:
    type: File[]
    inputBinding:
      position: 1
  min_length:
    type: int?
    inputBinding:
      prefix: "-l"
  run_check_m:
    type: boolean?
    inputBinding:
      prefix: "--run-checkm"
    default: false
  outdir:
    type: string?
    inputBinding:
      prefix: "-o"
    default: "output_metawrap"


outputs:
  concoct_bins:
    type: Directory?
    outputBinding:
      glob: $(inputs.outdir)/concoct_bins
  metabat2_bins:
    type: Directory?
    outputBinding:
      glob: $(inputs.outdir)/metabat2_bins
  maxbin2_bins:
    type: Directory?
    outputBinding:
      glob: $(inputs.outdir)/maxbin2_bins


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"

doc: |
  https://arxiv.org/abs/1604.03071
  http://cab.spbu.ru/files/release3.12.0/manual.html#meta