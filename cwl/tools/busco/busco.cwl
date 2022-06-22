cwlVersion: v1.2
class: CommandLineTool
label: BUSCO
doc: |
      Estimation of Eukaryotic MAG completeness

requirements:
  ResourceRequirement:
    coresMin: 4
    ramMin: 5000
hints:
  #DockerRequirement:
  #  dockerPull: todo

baseCommand: [ 'busco' ]

arguments:
  - --offline
  - --auto-lineage-euk
  - valueFrom: $(runtime.cores)
    prefix: "-c"

#check is --download_path option required? Should be mounted in container
inputs:
  mag:
    type: File
    label: MAG fasta
    format: edam:format_1929 #FASTA
    inputBinding:
      position: 1
      prefix: -i
  mode:
    type: string
    label: busco mode: genome, proteins or transcriptome
    inputBinding:
      position: 2
      prefix: -m
      default: 'genome'

outputs:
  summary:
    type: File
    outputBinding:
      glob: "short_summary.specific*.out.txt"


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20