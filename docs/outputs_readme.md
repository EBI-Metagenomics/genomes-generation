# Genomes Generation Pipeline output structure

GGP generates two types of outputs: **basic** or **expanded** results ([exp]), which additionally include the intermediate files produced during the main processing steps.
Expanded results can be included using argument `--publish_all`.

## Top level structure
```commandline
├── input
├── qc            [exp]
├── binning       [exp]
├── eukaryotes
├── prokaryotes
├── upload        [optional]
└── pipeline_info
```

### input
`samplesheet.csv`: generated pipeline input samplesheet \
`runs_assemblies.tsv`: mapping file contains raw_reads_run_identifier and corresponding assembly_identifier

### qc [expanded output]
Quality control statistics. 

TODO: add fastQC
```commandline
qc
├── fastp
│   ├── *fastp.html
│   ├── *fastp.json
│   └── *fastp.log
```

### binning [expanded output]
Binning results for CONCOCT, MetaBAT2 and Maxbin2. MaxBin2 and MetaBAT2 produce folder with `descarded` bins by several criteria. CONCOCT and MetaBAT2 also contain `coverage` files.

```commandline
binning
├── concoct
│   ├── <bins>*.fa
│   └── coverage.tsv
├── maxbin2
│   ├── <bins>*.fa
│   └── discarded
├── metabat2
│   ├── <bins>*.fa
│   ├── depth.tsv
│   └── discarded
```

### eukaryotes 
```commandline
eukatyotes
├── bins
│   ├── <bins>.fa.gz
├── coverage
│   ├── <mag>_coverage.txt
│   └── aggregated_contigs2bins.txt
├── mags
│   ├── <mags>.fa.gz
├── refinement                        [exp] 
│   ├── eukcc
│   │   ├── merged_bins.csv
│   │   └── 
│   ├── binlinks
├── stats
│   ├── eukcc_final_qc.csv
│   └── combined_busco_eukcc.csv
└── taxonomy
    ├── all_bin2classification.txt
    └── human_readble.tsv
```

### prokaryotes 
```commandline
prokaryotes
├── bins
│   ├── <bins>.fa.gz
├── coverage
│   ├── <mag>_coverage.txt
│   └── aggregated_contigs2bins.txt
├── mags
│   ├── <mags>.fa.gz
├── refinement                        [exp] 
│   ├── 
├── rna  
├── stats
│   ├── 
│   └── 
└── taxonomy
    ├── 
    └── 
```

### upload [optional]
```
    ├── unclassified_genomes.txt
    ├── final_table_for_uploader.tsv
    ├── upload
    │  ├── create_manifests
    │  │   ├── ENA_backup.json
    │  │   ├── genome_samples.xml
    │  │   ├── manifests
    │  │      ├── MAG_IDENTIFIER.manifest
    │  │      ├── registered_MAGs.tsv
    │  │      └── submission.xml
    │  ├── ena_submission_summary.txt
    │  └── webin_cli
    │      ├── MAG_IDENTIFIER_webin-cli.report

```

### pipeline_info
```
└── pipeline_info/
    ├── structured_pipeline_logging_by_runs.txt
    ├── execution_report_*.html
    ├── pipeline_dag_*.html
    └── software_versions.yml
```
- `software_versions.yml` tools versions used in pipeline
