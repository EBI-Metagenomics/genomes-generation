# EukRecover
Pipeline to recover eukaryotic MAGs using CONCOCT, metaBAT2 and EukCC's merging algorythm.

Needs paired end shotgun metagenomic reads.

## Environment

Eukrecover requires an environment with snakemake and metaWRAP.

## Quickstart

Define your samples in the file `samples.csv`.
This file needs to have the columns project and run to identify each metagenome. 

This pipeline does not support co-binning, but feel free to change it. 

Clone this repro wherever you want to run the pipeline:
```
git clone https://github.com/openpaul/eukrecover/
```


You can then run the snakemake like so

```
snakemake --use-singularity
```

The pipeline used dockerhub to fetch all tools, so make sure you have singularity installed.



## Prepare databases
The pipeline will setup databases for you, but if you already have a EukCC or a BUSCO 5 database you can use them 
by specifying the location in the file `config/config.yaml`


## Output:
In the folder results you will find a folder `MAGs` which will contain a folder
`fa` containing the actual MAG fastas.
In addition you will find stats for each MAG in the table `QC.csv`.

This table contains the following columns:

name,eukcc_compl,eukcc_cont,BUSCO_C,BUSCO_M,BUSCO_D,BUSCO_F,BUSCO_tax,N50,bp



## Citation:

If you use this pipeline please make sure to cite all used software. 

For this please reffer to the used rules.
