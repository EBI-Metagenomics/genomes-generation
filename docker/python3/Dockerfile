FROM python:3.7.9-slim-buster

LABEL maintainer="Microbiome Informatics Team www.ebi.ac.uk/metagenomics"

RUN apt-get update && rm -rf /var/lib/apt/lists/*

RUN /usr/local/bin/python -m pip install --upgrade pip && pip install --no-cache-dir -q biopython

COPY requirements.txt /
RUN pip install -r /requirements.txt

COPY binlinks.py \
     create_qc_table.py \
     rename_multifasta_prefix.py \
/tools/

RUN chmod a+x /tools/*

ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/tools