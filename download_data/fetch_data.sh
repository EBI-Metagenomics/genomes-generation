#!/usr/bin/env bash

function FetchReads {
  echo 'Fetching reads...'
  . "/hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh"
  mitload fetchtool
  bsub -I -q production "bash fetch-reads-tool.sh -v -p $READS_ACC -d ${CATALOGUE_PATH}/Raw_reads/" --run-list ${CATALOGUE_PATH}/runs.tsv
}


function FetchAssemblies {
  echo 'Fetching assemblies...'
  . "/hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh"
  mitload fetchtool
  bsub -I -q production "bash fetch-assemblies-tool.sh -v -p $SAMPLE -d ${CATALOGUE_PATH}/Assemblies/"
}


function RunRenamingScript {
  echo 'Starting renaming script...'
  . /hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh
  mitload assembly_pipeline
  sleep 1
  mkdir -p ${CATALOGUE_PATH}/Uploaded_Assembly_IDs

  python3 ${REPO_PATH}/rename-erz.py \
      -d ${CATALOGUE_PATH}/Assemblies/${SAMPLE}/raw/ \
      -o ${CATALOGUE_PATH}/Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt

  cat ${CATALOGUE_PATH}/Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt | tr ',' '\t' > ${CATALOGUE_PATH}/rename.tsv
  cut -f1 ${CATALOGUE_PATH}/rename.tsv > ${CATALOGUE_PATH}/runs.tsv
  export CONVERT=${CATALOGUE_PATH}/Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt
}


function Rename {
  echo "Rename assemblies"
  while read line; do
    if [[ $line == [SED]RR* ]]
    then
      export OLD=$(echo $line | cut -d ',' -f1)
      export NEW=$(echo $line | cut -d ',' -f2)
      mv ${CATALOGUE_PATH}/Assemblies/${SAMPLE}/raw/${NEW}.fasta.gz ${CATALOGUE_PATH}/Assemblies/${SAMPLE}/raw/${OLD}.fasta.gz
    fi
  done < $CONVERT
}


function GenerateSamplesheet {
  echo "Generate samplesheet"
  python3 ${REPO_PATH}/generate_samplesheet.py \
  -a ${CATALOGUE_PATH}/Assemblies/${SAMPLE}/raw \
  -r ${CATALOGUE_PATH}/Raw_reads/${READS_ACC}/raw \
  -n ${CATALOGUE_PATH}/rename.tsv \
  -o ${CATALOGUE_PATH}/pipeline_input_samplesheet.csv
  echo "Done: ${CATALOGUE_PATH}/pipeline_input_samplesheet.csv"
}

SAMPLE=""
READS_ACC=""
CATALOGUE_PATH_INPUT=""
SKIP_FETCH="false" # By default fetch step included
SCRIPT_PATH="$(readlink -f $0)"
REPO_PATH="$(dirname $SCRIPT_PATH)"

while getopts 'a:r:c:f:' flag; do
    case "${flag}" in
        a) SAMPLE="$OPTARG" ;;
        r) READS_ACC="$OPTARG" ;;
        c) CATALOGUE_PATH_INPUT="$OPTARG" ;;
        f) SKIP_FETCH="$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

CATALOGUE_PATH=$(readlink -f "$CATALOGUE_PATH_INPUT")

if [[ -z $SAMPLE || -z $READS_ACC || -z $CATALOGUE_PATH || -z $REPO_PATH ]]; then
    echo "Missing required options"
    exit 1
fi

mkdir -p "$CATALOGUE_PATH"

if [[ $SKIP_FETCH == 'false' ]]; then
  FetchAssemblies
fi

RunRenamingScript
Rename

if [[ $SKIP_FETCH == 'false' ]]; then
  FetchReads
fi

GenerateSamplesheet