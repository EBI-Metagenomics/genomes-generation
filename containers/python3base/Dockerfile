FROM quay.io/biocontainers/python:3.9--1

LABEL maintainer="Microbiome Informatics Team www.ebi.ac.uk/metagenomics"
LABEL software="python"
LABEL software.version="3.9"
LABEL description="Python base image, includes the requests package"
LABEL website="https://github.com/EBI-Metagenomics/eukrecover"
LABEL license="Apache 2.0"

RUN /usr/local/bin/python -m pip install --upgrade pip && \
    pip install --no-cache-dir \
    pysam==0.21.0 \
    biopython==1.81
