#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

inputs:
  input_bam: File
  input_bins: Directory
  eukcc_db: Directory


outputs:

  linktable_file:
    type: File
    outputSource: linktable_bins/linktable
  eukcc_out:
    type: File
    outputSource: eukCC/eukcc_dir

steps:

  linktable_bins:
    run: ../../tools/linktable/linktable.cwl
    in:
      bam: input_bam
      bins: input_bins
    out:
      - linktable


  eukCC:
    run: ../../tools/eukcc2/eukcc_binmerging.cwl
    in:
      links: linktable_bins/linktable
      outdir:
        source: input_bins
        valueFrom: "$(self.basename)_eukcc"
      db: eukcc_db
      bins: input_bins
    out:
      - eukcc_dir