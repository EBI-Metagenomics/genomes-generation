cwlVersion: v1.2
class: CommandLineTool
label: drep
doc: |
      Pairwise comparison a list of genomes and dereplication using ANI thresholds

requirements:
  ResourceRequirement:
    coresMin: 2
    ramMin: 5000
  InlineJavascriptRequirement: {}

hints:
  DockerRequirement:
    dockerPull: "quay.io/microbiome-informatics/genomes-pipeline.drep:v2"

baseCommand: [ 'dRep', 'dereplicate' ]

arguments:
  - valueFrom: $(runtime.cores)
    prefix: "-p"
  #- prefix: '-g'
  #  valueFrom: $(inputs.genomes.listing)
  #  position: 3

#check is --download_path option required? Should be mounted in container
inputs:
  genomes:
    type: File[] #Directory
    inputBinding:
      position: 3
      prefix: -g
  pani:
    type: float?
    label: ANI threshold to form primary (MASH) clusters
    inputBinding:
      position: 3
      prefix: -pa
    default: 0.80
  sani:
    type: float?
    label: ANI threshold to form secondary (MASH) clusters
    inputBinding:
      position: 4
      prefix: -sa
    default: 0.99
  cov:
    type: float?
    label: Coverage threshold - minimum level of overlap between genomes
    inputBinding:
      position: 5
      prefix: -nc
    default: 0.40
  cov_method:
    type: string?
    label: Method to calculate coverage of an alignment - larger or total
    inputBinding:
      position: 6
      prefix: -cm
    default: "larger"
  quality:
    type: File
    label: Quality information on the genomes
    inputBinding:
      position: 7
      prefix: --genomeInfo
  completeness:
    type: int?
    label: Minumum genome completeness
    inputBinding:
      position: 8
      prefix: -comp
    default: 49
  contamination:
    type: int?
    label: Maximum genome contamination
    inputBinding:
      position: 9
      prefix: -con
    default: 21
  drep_outfolder:
    type: string?
    inputBinding:
      position: 11
    default: "drep_out"


outputs: []
  #drep_genomes:
  #  type: Directory
  #  outputBinding:
   #   glob: "dereplicated_genomes"


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-29