#!/bin/bash
#SBATCH --job-name=bin_upload_ggp
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=bin_upload.out
#SBATCH --error=bin_upload.err 


export NFX_OPTS="-Xms=512m -Xmx=8g"
export NXF_ANSI_LOG=false

export NXF_TEMP="/hps/scratch/rdf/metagenomics/nf-scratch/"
export TMPDIR="/hps/scratch/rdf/metagenomics/nf-scratch/"
export NXF_SINGULARITY_CACHEDIR="/hps/nobackup/rdf/metagenomics/service-team/singularity-cache/"
export NXF_ANSI_LOG=false

nextflow run bin_upload_tsv.nf \
-profile slurm,ebi \
-entry BIN_UPLOAD_TSV
