## Run GGP using ENA data

ENA is very useful resource for MAGs/bins generation because it [structured](https://www.ebi.ac.uk/training/online/courses/ena-quick-tour/what-is-ena/how-is-ena-structured/) raw-reads and retrieved assemblies in study format.

For example, 
`raw-reads study` [SRP080008](https://www.ebi.ac.uk/ena/browser/view/PRJNA330077) was assembled and submitted as `assembly study` [ERP108081](https://www.ebi.ac.uk/ena/browser/view/PRJEB26108).

It is possible to submit generated MAGs to ENA **if you own the data**. That pipeline is automated to register and submit MAGs using [genome_uploader](https://github.com/EBI-Metagenomics/genome_uploader) to ENA mags layer under `assembly study` project. 
> [!NOTE]
> If you do not have full access to ENA studies - do not submit MAGs 

## Simple run without MAGs upload
```bash
nextflow run ebi-metagenomics/genomes-generation \
-profile `specify profile(s)` \
--ena_assembly_study_accession `assembly study` \
--ena_raw_reads_study_accession `raw-reads study` \
--outdir `full path to output directory`
```

## Run with MAGs upload
```bash
nextflow run ebi-metagenomics/genomes-generation \
-profile `specify profile(s)` \
--ena_assembly_study_accession `assembly study` \
--ena_raw_reads_study_accession `raw-reads study` \
--metagenome "chosen metagenome" \
--biomes "chosen biome,chosen feature,chosen material" \ 
--centre_name "name" \
--upload_tpa true \
--upload_mags true \
--upload_bins false \
--outdir `full path to output directory`
```

## Run with filtering input data
```bash
nextflow run ebi-metagenomics/genomes-generation \
-profile `specify profile(s)` \
--ena_assembly_study_accession `assembly study` \
--ena_raw_reads_study_accession `raw-reads study` \
--filter_samplesheet_by_scientific_name "marine sediment metagenome,sediment metagenome" \
--filter_samplesheet_by_environment_biome "cold marine sediment biome,warm marine sediment biome" \
--keep_metat_runs_in_samplesheet false \
--outdir `full path to output directory`
```

## Optional arguments

### input samplesheet generation
**By default GGP filters out AMPLICON and metaT runs.**
- `--filter_samplesheet_by_scientific_name (default=false)`: filter by **scientific_name** (like in ENA records) comma separated list
- `--filter_samplesheet_by_environment_biome (default=false)`: filter by **environment_biome** (like in ENA records) comma separated list
- `--keep_metat_runs_in_samplesheet (default=false)`: allow **METATRANSCRIPTOMIC** runs (use that argument if you are 100% sure that ENA has mislabeled runs `library_source`)

### ENA upload
#### Metagenome
- `--metagenome ` \
Manually choose the most appropriate metagenome from https://www.ebi.ac.uk/ena/browser/view/408169?show=tax-tree. \
For example, `marine metagenome`
#### Environment information
- `--biomes ` \
Comma-separated environment parameters in format: 
`"environment_biome,environment_feature,environment_material"` \
For example, `marine sediments,subtropical gyre,sinking marine particle`
#### genome_uploader
- `--centre_name`: name of the centre generating and uploading genomes
- `--upload_tpa`: if uploading TPA (Third PArty) generated genomes
- `--upload_force`: forces reset of submission xml and ENA backup
- `--upload_mags`: upload MAGs (Metagenome-Assembled Genomes)
- `--upload_bins`: upload bins
- `--test_upload`: use **test** server for upload, not live mode
