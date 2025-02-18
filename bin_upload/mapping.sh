#!/bin/bash

ref_fasta=$1

echo " ---> indexing assembly"
run_id=$(basename "$ref_fasta" | cut -d'_' -f1)
bin_id=$(basename "$ref_fasta" | cut -d'.' -f1)

read_fwd="$run_id"_1.fastq.gz
read_rev="$run_id"_2.fastq.gz

singularity run "$SINGULARITY_CACHEDIR"/quay.io-microbiome-informatics-bwa_metabat_concoct-2.2.1_2.16_1.1.0.img bwa-mem2 index "$ref_fasta"

echo " ---> mapping files to assembly"
singularity run "$SINGULARITY_CACHEDIR"/quay.io-microbiome-informatics-bwa_metabat_concoct-2.2.1_2.16_1.1.0.img bwa-mem2 mem -M \
    -t 6 \
    "$ref_fasta" \
    "$read_fwd" "$read_rev" | \
singularity run "$SINGULARITY_CACHEDIR"/quay.io-microbiome-informatics-bwa_metabat_concoct-2.2.1_2.16_1.1.0.img samtools view -@ 6 -q 20 -Sb -F 4 - | \
singularity run "$SINGULARITY_CACHEDIR"/quay.io-microbiome-informatics-bwa_metabat_concoct-2.2.1_2.16_1.1.0.img samtools sort -@ 6 -O bam - -o "$bin_id"_sorted.bam
echo " ---> samtools index sorted bam"
singularity run "$SINGULARITY_CACHEDIR"/quay.io-microbiome-informatics-bwa_metabat_concoct-2.2.1_2.16_1.1.0.img samtools index -@ 6 "$bin_id"_sorted.bam

echo " ---> samtools idxstats sorted bam"
singularity run "$SINGULARITY_CACHEDIR"/quay.io-microbiome-informatics-bwa_metabat_concoct-2.2.1_2.16_1.1.0.img samtools idxstats --threads 6 "$bin_id"_sorted.bam > "$bin_id".assembly.idxstats

echo " ---> depth generation"
singularity run "$SINGULARITY_CACHEDIR"/quay.io-microbiome-informatics-bwa_metabat_concoct-2.2.1_2.16_1.1.0.img jgi_summarize_bam_contig_depths \
    --outputDepth "$bin_id".txt \
    "$bin_id"_sorted.bam

echo " --> removing tmp files"
