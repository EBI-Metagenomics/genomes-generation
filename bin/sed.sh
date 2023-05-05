#!/bin/bash
while getopts i: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
    esac
done

echo "Change contig headers"
sed 's/\./_/' ${INPUT} > out
echo "Done"