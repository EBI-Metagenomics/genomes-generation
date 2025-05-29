#!/usr/bin/env bash

export SCIENTIFIC_NAME=""
export FILTER_METAT=""

while getopts 'a:r:c:b:f:' flag; do
    case "${flag}" in
        a) ASSEMBLY_STUDY_ACC="$OPTARG" ;;
        r) RAW_READS_STUDY_ACC="$OPTARG" ;;
        c) CATALOGUE_PATH_INPUT="$OPTARG" ;;
        b) SCIENTIFIC_NAME="$OPTARG" ;;
        f) FILTER_METAT="$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

export SCRIPT_PATH="$(readlink -f $0)"
export REPO_PATH="$(dirname $SCRIPT_PATH)"

. /hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh
mitload assembly_pipeline
sleep 1

args=(
  -a "$ASSEMBLY_STUDY_ACC"
  -r "$RAW_READS_STUDY_ACC"
  -o "$CATALOGUE_PATH_INPUT"
)

if [ -n "$SCIENTIFIC_NAME" ]; then
  echo "Scientific name: $SCIENTIFIC_NAME"
  args+=(--scientific-name "$SCIENTIFIC_NAME")
fi

if [ -n "$FILTER_METAT" ]; then
  echo "Filtering is on"
  args+=(--filter-out-metat)
fi

echo "Generate assemblies and runs filtered list"

python3 ${REPO_PATH}/filter_runs.py "${args[@]}"

cut -f1 "${CATALOGUE_PATH_INPUT}/runs_assemblies.tsv" > ${CATALOGUE_PATH_INPUT}/runs.tsv
cut -f2 "${CATALOGUE_PATH_INPUT}/runs_assemblies.tsv" > ${CATALOGUE_PATH_INPUT}/assemblies.tsv

echo "Done: ${CATALOGUE_PATH_INPUT}/runs_assemblies.tsv"
