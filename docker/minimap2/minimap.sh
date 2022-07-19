#!/bin/bash

set -e

usage()
{
cat << EOF
usage: $0 options

Map cleaned reads to assembly

OPTIONS:
   -s      Scaffolds fasta file
   -f      Forward or single-end fastq file (.fastq or .fastq.gz) [REQUIRED]
   -r	     Reverse fastq file (.fastq or .fastq.gz) [OPTIONAL]
   -t      Number of threads [OPTIONAL]
EOF
}

while getopts :s:f:r:t: option; do
    case "${option}" in
        s) SCAFFOLDS=${OPTARG};;
        f) FASTQ_R1=${OPTARG};;
        r) FASTQ_R2=${OPTARG};;
        t) THREADS=${OPTARG};;
        *) echo "invalid option"; exit;;
    esac
done

# check if all required arguments are supplied
if [[ -z ${FASTQ_R1} ]]
then
     echo "ERROR : Please supply input FASTQ files"
     usage
     exit 1
fi

if [ ${THREADS} -eq 1 ]
then
    THREADS_SAM=1
else
    THREADS_SAM=$((${THREADS}-1))
fi


# ?? option for single end ??

filename_scaffolds=$(basename -- "$SCAFFOLDS")
name_scaffolds=${filename_scaffolds%_*}

echo "mapping reads to contigs ${name_scaffolds}"
minimap2 -ax sr -t $THREADS_SAM $SCAFFOLDS $FASTQ_R1 $FASTQ_R2 |
samtools view -q 20 -Sb - | \
samtools sort -@ $THREADS_SAM -O bam - -o ${name_scaffolds}.bam
echo "index reads"
samtools index ${name_scaffolds}.bam

#rm ${filename_scaffolds} ${filename_reads1} ${filename_reads2}








