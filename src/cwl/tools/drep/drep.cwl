cwlVersion: v1.2
class: CommandLineTool
label: drep
doc: |
      Pairwise comparison a list of genomes and dereplication using ANI thresholds

requirements:
  ResourceRequirement:
    coresMin: 4
    ramMin: 5000
hints:
  #DockerRequirement:
  #  dockerPull: todo

baseCommand: [ 'dRep', 'dereplicate', 'working_directory' ]

arguments:
  - valueFrom: $(runtime.cores)
    prefix: "-p"

#check is --download_path option required? Should be mounted in container
inputs:
  genomes:
    type: File
    label: text file with paths of genomes
    format: edam:format_1964
    inputBinding:
      position: 1
      prefix: -g
  pani:
    type: float
    label: ANI threshold to form primary (MASH) clusters
    inputBinding:
      position: 2
      prefix: -pa
      default: 0.80
  sani:
    type: float
    label: ANI threshold to form secondary (MASH) clusters
    inputBinding:
      position: 3
      prefix: -sa
      default: 0.95
  cov:
    type: float
    label: Coverage threshold: minimum level of overlap between genomes
    inputBinding:
      position: 4
      prefix: -nc
      default: 0.40
  cov_method:
    type: string
    label: Method to calculate coverage of an alignment: larger or total
    inputBinding:
      position: 5
      prefix: -cm
      default: "larger"
  quality:
    type: File
    label: Quality information on the genomes
    format: format_3752
    inputBinding:
      position: 6
      prefix: --genomeInfo
  completeness:
    type: int
    label: Minumum genome completeness
    inputBinding:
      position: 7
      prefix: -comp
      default: 49
  contamination:
    type: int
    label: Maximum genome contamination
    inputBinding:
      position: 8
      prefix: -con
      default: 21

outputs:
  drep_genomes:
    type: Directory
    outputBinding:
      glob: working_directory


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20