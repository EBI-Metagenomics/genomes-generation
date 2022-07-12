#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}


inputs:
  concoct_csv: File
  metabat2_csv: File
  concoct_bins: Directory
  metabat2_bins: Directory

outputs:
  dereplicated_genomes:
    type: Directory
    outputSource: dRep/drep_genomes

steps:

  choose_genomes:
    run: ../../tools/drep/choose_genomes.cwl
    in:
      concoct_csv: concoct_csv
      metabat2_csv: metabat2_csv
      concoct_merged: concoct_bins
      metabat2_merged: metabat2_bins
    out:
      - quality
      - input_genomes

  dRep:
    run: ../../tools/drep/drep.cwl
    in:
      genomes: choose_genomes/input_genomes
      quality: choose_genomes/quality
    out:
      - drep_genomes


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20