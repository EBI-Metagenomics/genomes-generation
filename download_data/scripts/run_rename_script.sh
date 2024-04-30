#!/usr/bin/env bash

while getopts 'a:c:' flag; do
    case "${flag}" in
        a) ASSEMBLY_STUDY_ACC="$OPTARG" ;;
        c) CATALOGUE_PATH_INPUT="$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

export SCRIPT_PATH="$(readlink -f $0)"
export REPO_PATH="$(dirname $SCRIPT_PATH)"

. /hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh
mitload assembly_pipeline
sleep 1

echo "Generate rename"

python3 ${REPO_PATH}/rename-erz.py \
  -d ${CATALOGUE_PATH_INPUT}/Assemblies/${ASSEMBLY_STUDY_ACC}/raw \
  -o ${CATALOGUE_PATH_INPUT}/rename.tsv

cut -f1 ${CATALOGUE_PATH_INPUT}/rename.tsv > ${CATALOGUE_PATH_INPUT}/runs.tsv

echo "Done: ${CATALOGUE_PATH_INPUT}/rename.tsv"
