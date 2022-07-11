cwlVersion: v1.2
class: CommandLineTool
label: EukCC2 bin merging
doc: |
      Identify medium quality bins and merge based on linked reads and improved completeness/contamination thresholds

requirements:
  #InitialWorkDirRequirement:
  #  listing:
  #    - $(inputs.bins)
  ResourceRequirement:
    coresMin: 4
    ramMin: 5000
hints:
  DockerRequirement:
    dockerPull: "quay.io/microbiome-informatics/eukcc:latest"

baseCommand: [ 'eukcc' ]

arguments:
  - valueFrom: debug_$(inputs.bins.basename)
    prefix: "--debug"
  - valueFrom: ".fa"
    prefix: "--suffix"
  - valueFrom: $(runtime.cores)
    prefix: "--threads"
  - valueFrom: $(inputs.links.basename)_merged.
    prefix: "--prefix"

#check is --db option required? Should be mounted in container
inputs:
  improve_percent:
    type: int?
    label: 100-n% completeness to constitute primary bin
    inputBinding:
      position: 1
      prefix: --improve_percent
    default: 10
  n_combine:
    type: int?
    label: compare n primary bins with secondary bins
    inputBinding:
      position: 2
      prefix: --n_combine
    default: 1
  improve_ratio:
    type: int?
    label: increased n * contamination should be less than gained completeness of merged bin
    inputBinding:
      position: 3
      prefix: --improve_ratio
    default: 5
  min_links:
    type: int?
    label: find bins linked by at least n number of paired reads
    inputBinding:
      position: 4
      prefix: --min_links
    default: 100
  links:
    type: File
    label: csv of linking reads
    format: edam:format_3752
    inputBinding:
      position: 5
      prefix: --links
  outdir:
    type: string
    label: output directory path
    inputBinding:
      position: 6
      prefix: --out
  db:
    type: Directory
    label: directory containing eukcc database
    inputBinding:
      position: 7
      prefix: --db
  bins:
    type: Directory
    label: directory containing bins with ext *.fa
    inputBinding:
      position: 8

outputs:
  eukcc_dir:
    type: Directory
    outputBinding:
      glob: "$(inputs.outdir)"
  csv:
    type: File
    outputBinding:
      glob: "$(inputs.outdir)/eukcc.csv"


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
s:dateCreated: 2022-06-20