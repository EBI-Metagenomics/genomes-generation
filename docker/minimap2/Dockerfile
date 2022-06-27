FROM alpine:3.7 as build

LABEL maintainer="Microbiome Informatics Team www.ebi.ac.uk/metagenomics"
LABEL software="minimap2"
LABEL description="Map reads to scaffolds and select high quality alignments. Sort and output bam file."

RUN apk add --no-cache build-base ncurses-dev bzip2-dev xz-dev zlib-dev

ENV SAMTOOLS_VERSION=1.9
# samtools and bgzip
RUN mkdir /samtools && \
    wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 && \
    tar -xjf samtools-$SAMTOOLS_VERSION.tar.bz2 && rm samtools-$SAMTOOLS_VERSION.tar.bz2 && \
    cd samtools-$SAMTOOLS_VERSION && ./configure --prefix=/samtools && make && make install install-htslib

ENV MINIMAP2_VERSION=2.24
RUN wget https://github.com/lh3/minimap2/releases/download/v$MINIMAP2_VERSION/minimap2-$MINIMAP2_VERSION.tar.bz2 && \
    tar -xjf minimap2-$MINIMAP2_VERSION.tar.bz2 && rm minimap2-$MINIMAP2_VERSION.tar.bz2 && \
    cd minimap2-$MINIMAP2_VERSION && make


FROM alpine:3.7
ENV MINIMAP2_VERSION=2.24
RUN apk add --no-cache pigz nodejs coreutils bash zlib bzip2 ncurses  \
    ncurses-dev bzip2-dev xz-dev zlib-dev tar && \
    mkdir /tools /samtools /minimap2

COPY --from=build /samtools /samtools
COPY --from=build minimap2-$MINIMAP2_VERSION /minimap2

COPY minimap.sh /tools/

RUN chmod a+x /tools/*

ENV PATH="/samtools/bin:/tools:/minimap2:${PATH}"