#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "metaWrap binning tool (binrefine)"

hints:
  DockerRequirement:
    dockerPull: "quay.io/microbiome-informatics/metawrap:latest"

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 8
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.refdata)
        entryname: $("/dbs/checkm_db")

baseCommand: ['metawrap', 'bin_refinement', '--quick']

arguments:
  - valueFrom: $(runtime.ram)
    prefix: -m
  - valueFrom: $(runtime.cores)
    prefix: -t

inputs:
  concoct_bin_dir:
    type: Directory?
    inputBinding:
      prefix: "-A"
  metabat_bin_dir:
    type: Directory?
    inputBinding:
      prefix: "-B"
  maxbin_bin_dir:
    type: Directory?
    inputBinding:
      prefix: "-C"
  completion:
    type: int
    default: 50
    inputBinding:
      prefix: "-c"
  contamination:
    type: int
    default: 10
    inputBinding:
      prefix: "-x"
  out_dir:
    type: string?
    default: "metaWrap_binrefine"
    inputBinding:
      prefix: "-o"
  refdata:
    type: Directory


outputs:
  bins:
    type: Directory
    outputBinding:
      glob: $(inputs.out_dir)


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