cwlVersion: v1.2
class: CommandLineTool
label: Create table of linking reads
doc: |
      Select reads aligning to bin with criteria: minimum ANI & contig mapping location

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.bins)
  ResourceRequirement:
    coresMin: 4
    ramMin: 5000
hints:
  #DockerRequirement:
  #  dockerPull: todo

baseCommand: [ 'python', 'binlinks.py' ]

arguments:
  - valueFrom: $(inputs.bam.basename).links.csv
    prefix: "--out"

inputs:
  ani:
    type: int?
    label: minimum average nucletoide identity
    inputBinding:
      position: 1
      prefix: --ANI
      default: 99
  within:
    type: int
    label: keep reads mapping within the first or last n basepairs of contig
    inputBinding:
      position: 2
      prefix: --within
      default: 1500
  bam:
    type: string
    label: output filename
    inputBinding:
      position: 3
  bins:
    type: Directory
    label: directory containing bins with ext *.fa
    inputBinding:
      position: 4

#edit if outreads1 there should be outreads2
outputs:
  linktable:
    type: File
    format: edam:format_3752
    outputBinding:
      glob: "*.links.csv"


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20