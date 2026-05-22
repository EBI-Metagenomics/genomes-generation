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
`assembly_software.tsv`: file listing assembler tool used to assemble each run \
`runs_assemblies.tsv`: mapping file contains raw_reads_run_identifier and corresponding assembly_identifier

### qc [expanded output]
Quality control statistics. 

TODO: add fastQC
```commandline
qc
└── fastp
    └── <run>
        ├── <run>.fastp.html
        ├── <run>.fastp.json
        └── <run>.fastp.log
```

### binning [expanded output]
Binning results for CONCOCT, MetaBAT2 and Maxbin2. MaxBin2 and MetaBAT2 produce folder with `descarded` bins by several criteria. CONCOCT and MetaBAT2 also contain `coverage` files.

```commandline
binning
├── concoct
│   └── <run>
│       ├── bins
│       │   └── <bin>.fa
│       └── coverage
│           └── *.tsv
├── maxbin2
│   └── <run>
│       ├── bins
│       │   └── <bin>.fa
│       └── discarded
└── metabat2
    └── <run>
        ├── bins
        │   └── <bin>.fa
        ├── depth
        │   └── *.txt.gz
        └── discarded

### eukaryotes 
```commandline
eukatyotes
├── bins                                                   [exp]
│   └── <bin>.fa.gz
├── coverage
│   ├── <mag>_coverage.txt
│   └── <project>_contigs2bins.txt
├── drep
│   ├── per-run                                            [exp]
│   │   └── <run>
│   │       ├── data_tables
│   │       │   └── *.csv
│   │       └── dereplicated_genomes.tsv
│   ├── data_tables                                        [exp] 
│   │   └── *.csv
│   └── <project>_dereplicated_genomes.tsv
├── mags
│   └── <mag>.fa.gz
├── refinement                                             [exp] 
│   ├── eukcc
│   │   └── <run>
│   │       ├── merged_bins.csv
│   │       └── eukcc.csv
│   └── binlinks
│       └── <run>.links.csv
├── stats
│   ├── <project>_eukcc_before_filter_and_dedup.csv
│   ├── <project>_eukcc_bins_quality_filtered.csv
│   └── <project>_eukcc_busco_bins_quality_filtered.csv      [exp] 
└── taxonomy
    ├── <project>_bins_ncbi_taxonomy.txt
    ├── <project>_mags_bat_output.txt
    └── <project>_mags_ncbi_taxonomy.txt
```

### prokaryotes 
```commandline
prokaryotes
├── bins                                                   [exp]
│   └── <bin>.fa.gz
├── coverage
│   ├── <mag>_coverage.txt
│   └── <project>_contigs2bins.txt
├── drep
│   ├── data_tables                                        [exp] 
│   │   └── *.csv
│   └── <project>_dereplicated_genomes.tsv
├── mags
│   └── <mag>.fa.gz
├── refinement                                             [exp] 
│   └── binette
│       ├── <run>_final_bins_quality_reports.tsv
│       └── <run>_input_bins_quality_reports
│           └── <binner>*.tsv
├── rna
│   └── <bin>
│       ├── <bin>_rRNAs.fasta
│       ├── <bin>_rRNAs.out
│       └── <bin>_tRNA_20aa.out
├── stats
│   ├── <project>_checkm2_all_bins.csv
│   ├── <project>_checkm2_bins_quality_filtered.tsv
│   └── <project>_gunc_contamination_report.txt
└── taxonomy
    ├── <project>_bins_ncbi_taxonomy.txt
    ├── <project>_gtdbtk_results.tar.gz
    └── <project>_mags_ncbi_taxonomy.txt
```

### upload [optional output]
```
└── upload
    ├── bins
    │   ├── ena_submission_summary.txt
    │   ├── genome_uploader
    │   │   ├── ENA_backup.json
    │   │   ├── final_table_for_uploader.tsv
    │   │   ├── genome_samples.xml
    │   │   ├── manifests
    │   │   │   └── <bin>.manifest
    │   │   └── registered_bins.tsv
    │   └── webin_cli
    │       └── *_webin-cli.report
    └── mags
        ├── ena_submission_summary.txt
        ├── genome_uploader
        │   ├── ENA_backup.json
        │   ├── final_table_for_uploader.tsv
        │   ├── genome_samples.xml
        │   ├── manifests
        │   │   └── <mag>.manifest
        │   └── registered_MAGs.tsv
        └── webin_cli
            └── *_webin-cli.report
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
