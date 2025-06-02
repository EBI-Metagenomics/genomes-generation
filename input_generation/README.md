## Pipeline input files generation for ENA data

If you are using ENA studies - use [generation script](input_generation/generate_inputs.py) to get both required input files. You need to provide `assembly_study_accession` and `raw_reads_study_accession`.

**! By default GGP filters out AMPLICON runs and metaT runs.**

Additional filters:
- filter by scientific_name use `-b SCIENTIFIC_NAME` (like in ENA records) comma separated list
- allow METATRANSCRIPTOMIC runs with `--keep-metat` (use that argument if you are 100% sure that ENA has mislabeled runs `library_source`)
```commandline
python3 input_generation/generate_inputs.py \
    -a assembly_study_accession \
    -r raw_reads_study_accession \
    [ -o output_directory_path ] \
    [ --output-samplesheet input_samplesheet.csv ] \
    [ --assembly-software-filename input_software.tsv ] \
    [ -b "marine sediment metagenome,sediment metagenome" ] \
    [ --keep-metat ]
```

*Do not forget to identify `metagenome` and `environment information` for pipeline execution*