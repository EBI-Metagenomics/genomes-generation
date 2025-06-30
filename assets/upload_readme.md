## Inputs for upload step

There are additional pipeline input arguments:

### assembly_software.tsv

`id`: samplesheet identifier \
`assembly_software`: tool that was used to assemble reads into contigs.

| id  | assembly_software  |
|-----|--------------------|
| ID  | assembler_vVersion |

> [!NOTE]
> If you used `generate_inputs.py` from [instuctions](../input_generation/README.md) then `assembly_software.tsv` should be already generated.

### Metagenome
Manually choose the most appropriate metagenome from https://www.ebi.ac.uk/ena/browser/view/408169?show=tax-tree. \
For example, `marine metagenome`

### Environment information
Comma-separated environment parameters in format: 
`"environment_biome,environment_feature,environment_material"` \
For example, `marine sediments,subtropical gyre,sinking marine particle`

## Run pipeline

```bash
nextflow run ebi-metagenomics/genomes-generation \
-profile `specify profile(s)` \
--input `samplesheet.csv` \
--assembly_software_file `software.tsv` \
--metagenome "chosen metagenome" \
--biomes "chosen biome,chosen feature,chosen material" \
--outdir `full path to output directory`
```

