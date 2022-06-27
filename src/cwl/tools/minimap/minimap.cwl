cwlVersion: v1.2
class: CommandLineTool
label: Mapping with minimap and samtools
doc: |
      Map reads to scaffolds and select high quality alignments. Sort and output bam file.

requirements:
  ResourceRequirement:
    coresMin: 4
    ramMin: 5000
hints:
  #DockerRequirement:
  #  dockerPull: todo

baseCommand: [ 'minimap.sh' ]

arguments:
- -t
- $(runtime.cores)

inputs:
  reads1:
    type: File
    format: edam:format_1930  # FASTQ
    label: cleaned reads forward file
    inputBinding:
      position: 1
      prefix: -f
  reads2:
    type: File?
    format: edam:format_1930  # FASTQ
    label: cleaned reads reverse file
    inputBinding:
      position: 2
      prefix: -r
  scaffolds:
    type: File
    format: edam:format_1929 #FASTA
    inputBinding:
      position: 3
      prefix: -s

#edit if outreads1 there should be outreads2
outputs:
  bamfile:
    type: File
    format: edam:format_2572
    secondaryFiles: .bai
    outputBinding:
      glob: "*.bam"


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20