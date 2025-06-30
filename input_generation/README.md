## Pipeline input files generation for ENA data

ENA is very useful resource for MAGs/bins generation because it [structured](https://www.ebi.ac.uk/training/online/courses/ena-quick-tour/what-is-ena/how-is-ena-structured/) raw-reads and retrieved assemblies in study format.

For example, 
raw-reads study [SRP080008](https://www.ebi.ac.uk/ena/browser/view/PRJNA330077) was assembled and submitted as [ERP108081](https://www.ebi.ac.uk/ena/browser/view/PRJEB26108).

If you are using ENA studies - use [generation script](input_generation/generate_inputs.py) to get both required input files. You need to provide `assembly_study_accession` and `raw_reads_study_accession`.

**! By default GGP filters out AMPLICON and metaT runs.**

Additional filters:
- filter by **scientific_name** use `-b SCIENTIFIC_NAME` (like in ENA records) comma separated list
- filter by **environment_biome** use `-e ENVIROMENT_BIOME` (like in ENA records) comma separated list
- allow **METATRANSCRIPTOMIC** runs with `--keep-metat` (use that argument if you are 100% sure that ENA has mislabeled runs `library_source`)
```commandline
python3 input_generation/generate_inputs.py \
    -a assembly_study_accession \
    -r raw_reads_study_accession \
    [ -o output_directory_path ] \
    [ --output-samplesheet input_samplesheet.csv ] \
    [ --assembly-software-filename input_software.tsv ] \
    [ -b "marine sediment metagenome,sediment metagenome" ] \
    [ -e "cold marine sediment biome,warm marine sediment biome" ] \
    [ --keep-metat ]
```