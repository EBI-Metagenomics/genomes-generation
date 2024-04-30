#!/usr/bin/env bash

while getopts 'a:c:r:' flag; do
    case "${flag}" in
        a) ASSEMBLY_STUDY_ACC="$OPTARG" ;;
        c) CATALOGUE_PATH="$OPTARG" ;;
        r) READS_STUDY_ACC="$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

export SCRIPT_PATH="$(readlink -f $0)"
export REPO_PATH="$(dirname $SCRIPT_PATH)"

. /hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh
mitload assembly_pipeline
sleep 1

echo "Generate samplesheet"

python3 ${REPO_PATH}/generate_samplesheet.py \
  -a ${CATALOGUE_PATH}/Assemblies/${ASSEMBLY_STUDY_ACC}/raw \
  -r ${CATALOGUE_PATH}/Raw_reads/${READS_STUDY_ACC}/raw \
  --assembly-software-filename ${CATALOGUE_PATH}/per_run_assembly.tsv \
  -o ${CATALOGUE_PATH}/pipeline_input_samplesheet.csv

echo "Done: ${CATALOGUE_PATH}/pipeline_input_samplesheet.csv"
