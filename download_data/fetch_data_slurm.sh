#!/usr/bin/env bash

function FilterRuns {
  echo 'Getting runs and assemblies'
  export filter_job_id=$(sbatch \
  -J filter_$ASSEMBLY_STUDY_ACC \
  --time=1:00:00 \
  --mem=5G \
  --mail-user=$USER@ebi.ac.uk \
  --mail-type=ALL \
  --wrap="bash ${REPO_PATH}/scripts/run_filter_script.sh -a $ASSEMBLY_STUDY_ACC -r $READS_STUDY_ACC -c $CATALOGUE_PATH -b ${SCIENTIFIC_NAME} -f ${FILTER_METAT}" | awk '{print $NF}')
  echo "Submitted FilterRuns $filter_job_id"
}

function FetchAssembliesSlurm {
  echo 'Fetching assemblies on Slurm...'
  export fetch_assemblies_job_id=$(sbatch --dependency=afterok:${filter_job_id} -J fetch_assemblies_$ASSEMBLY_STUDY_ACC --time=5-00:00:00 --mem=5G --mail-user=$USER@ebi.ac.uk --mail-type=ALL --wrap="bash fetch-assemblies-tool.sh -v -p $ASSEMBLY_STUDY_ACC -d ${CATALOGUE_PATH}/Assemblies/ --assembly-list ${CATALOGUE_PATH}/assemblies.tsv " | awk '{print $NF}')
  echo "Submitted FetchAssembliesSlurm $fetch_assemblies_job_id"
}

function FetchReadsSlurm {
  echo 'Fetching reads on Slurm...'
  export fetch_reads_job_id=$(sbatch --dependency=afterok:${filter_job_id} -J fetch_reads_$READS_STUDY_ACC --time=5-00:00:00 --mem=5G --mail-user=$USER@ebi.ac.uk --mail-type=ALL --wrap="bash fetch-reads-tool.sh -v -p $READS_STUDY_ACC -d ${CATALOGUE_PATH}/Raw_reads/ --run-list ${CATALOGUE_PATH}/runs.tsv" | awk '{print $NF}')
  echo "Submitted FetchReadsSlurm $fetch_reads_job_id"
}

export dependency_condition=""
function GenerateSamplesheet {
  if [ -n "$fetch_reads_job_id" ]; then
    export dependency_condition="--dependency=afterok:${fetch_reads_job_id}"
  fi

  if [ -n "$fetch_assemblies_job_id" ]; then
    if [ -n "$dependency_condition" ]; then
      export dependency_condition="${dependency_condition}:${fetch_assemblies_job_id}"
    else
      export dependency_condition="--dependency=afterok:${fetch_assemblies_job_id}"
    fi
  fi

  sbatch ${dependency_condition} \
  -J samplesheet_$ASSEMBLY_STUDY_ACC \
  --time=1:00:00 \
  --mem=5G \
  --mail-user=$USER@ebi.ac.uk \
  --mail-type=ALL \
  --wrap="bash ${REPO_PATH}/scripts/run_samplesheet_generation.sh -a $ASSEMBLY_STUDY_ACC -c $CATALOGUE_PATH -r $READS_STUDY_ACC"
}

export ASSEMBLY_STUDY_ACC=""
export READS_STUDY_ACC=""
export CATALOGUE_PATH_INPUT=""
export SCIENTIFIC_NAME=""
export FILTER_METAT="true"
export SKIP_FETCH="false" # By default fetch step included
export SCRIPT_PATH="$(readlink -f $0)"
export REPO_PATH="$(dirname $SCRIPT_PATH)"

while getopts 'a:r:c:f:b:t:' flag; do
    case "${flag}" in
        a) ASSEMBLY_STUDY_ACC="$OPTARG" ;;
        r) READS_STUDY_ACC="$OPTARG" ;;
        c) CATALOGUE_PATH_INPUT="$OPTARG" ;;
        f) SKIP_FETCH="$OPTARG" ;;
        b) SCIENTIFIC_NAME="$OPTARG" ;;
        t) FILTER_METAT="$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

export CATALOGUE_PATH=$(readlink -f "$CATALOGUE_PATH_INPUT")

if [[ -z $ASSEMBLY_STUDY_ACC || -z $READS_STUDY_ACC || -z $CATALOGUE_PATH || -z $REPO_PATH ]]; then
    echo "Missing required options"
    exit 1
fi


mkdir -p "$CATALOGUE_PATH"
if [[ $SKIP_FETCH == 'false' ]]; then
  . "/hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh"
  mitload fetchtool
  FilterRuns
  FetchAssembliesSlurm
  FetchReadsSlurm
fi

GenerateSamplesheet