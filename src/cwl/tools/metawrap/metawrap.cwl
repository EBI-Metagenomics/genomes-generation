#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: "metaWrap binning tool"

requirements:
#  DockerRequirement:
#    dockerPull: "quay.io/biocontainers/metawrap:1.1--0"
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: $(inputs.reads)
  ResourceRequirement:
    coresMin: 4

baseCommand: [ 'metawrap', 'binning', '--concoct', '--metabat2']

arguments:
  - valueFrom: $(runtime.outdir)
    prefix: -o
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


outputs:
  concoct_bins:
    type: Directory?
    outputBinding:
      glob: $("concoct_bins")
  metabat2_bins:
    type: Directory?
    outputBinding:
      glob: $("metabat2_bins")
  maxbin2_bins:
    type: Directory?
    outputBinding:
      glob: $("maxbin2_bins")
$namespaces:
 iana: https://www.iana.org/assignments/media-types/
 s: http://schema.org/

$schemas:
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"

doc: |
  https://arxiv.org/abs/1604.03071
  http://cab.spbu.ru/files/release3.12.0/manual.html#meta