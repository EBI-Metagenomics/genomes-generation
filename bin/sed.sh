#!/bin/bash
while getopts i: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
    esac
done

sed 's/\./_/' ${INPUT} > out