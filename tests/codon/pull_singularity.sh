#!/usr/bin/env bash

cd /hps/nobackup/rdf/metagenomics/singularity_cache


singularity pull --name quay.io_microbiome-informatics_eukcc:latest.sif \
docker://quay.io/microbiome-informatics/eukcc:latest

singularity pull --name quay.io_microbiome-informatics_eukrecover.python3_scripts:v1.sif \
docker://quay.io/microbiome-informatics/eukrecover.python3_scripts:v1

singularity pull --name quay.io_microbiome-informatics_metawrap:latest.sif \
docker://quay.io/microbiome-informatics/metawrap:latest

singularity pull --name quay.io_microbiome-informatics_eukrecover.minimap2.sh:v2.sif \
docker://quay.io/microbiome-informatics/eukrecover.minimap2.sh:v2

#&& mv pull/quay.io_microbiome-informatics_eukrecover.${NAME}:${VERSION}.sif .
