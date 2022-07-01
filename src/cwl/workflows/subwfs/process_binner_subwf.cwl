#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

doc: |
     Process binner results with linkmerge and eukCC2


requirements:
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  bam: File
  bins: Directory
  eukcc_db: Directory


outputs:
  linktable_file:
    type: File
    outputSource: linktable/linktable
  eukcc_out:
    type: File
    outputSource: eukCC/eukcc_dir


steps:

  linktable:
    run: ../../tools/linktable/linktable.cwl
    in:
      bam: bam
      bins: bins
    out:
      - linktable

  eukCC:
    run: ../../tools/eukcc2/eukcc_binmerging.cwl
    in:
      links: linktable/linktable
      outdir:
        source: bins
        valueFrom: $(self.basename)_eukcc
      db: eukcc_db
      bins: bins
    out:
      - eukcc_dir

$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20