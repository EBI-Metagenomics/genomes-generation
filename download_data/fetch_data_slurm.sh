#!/usr/bin/env bash

function FetchAssembliesSlurm {
  echo 'Fetching assemblies on Slurm...'
  export fetch_assemblies_job_id=$(sbatch -J fetch_assemblies_$ASSEMBLY_STUDY_ACC --time=5-00:00:00 --mem=5G --mail-user=$USER@ebi.ac.uk --mail-type=ALL --wrap="bash fetch-assemblies-tool.sh -v -p $ASSEMBLY_STUDY_ACC -d ${CATALOGUE_PATH}/Assemblies/" | awk '{print $NF}')
  echo "Submitted FetchAssembliesSlurm $fetch_assemblies_job_id"
}

function GenerateRename {
  export rename_job_id=$(sbatch --dependency=afterok:${fetch_assemblies_job_id} \
  -J rename_$ASSEMBLY_STUDY_ACC \
  --time=1:00:00 \
  --mem=5G \
  --mail-user=$USER@ebi.ac.uk \
  --mail-type=ALL \
  --wrap="bash ${REPO_PATH}/scripts/run_rename_script.sh -a $ASSEMBLY_STUDY_ACC -c $CATALOGUE_PATH" | awk '{print $NF}')
  echo "Submitted GenerateRename $rename_job_id"
}

function FetchReadsSlurm {
  echo 'Fetching reads on Slurm...'
  export fetch_reads_job_id=$(sbatch --dependency=afterok:${rename_job_id} -J fetch_reads_$READS_STUDY_ACC --time=5-00:00:00 --mem=5G --mail-user=$USER@ebi.ac.uk --mail-type=ALL --wrap="bash fetch-reads-tool.sh -v -p $READS_STUDY_ACC -d ${CATALOGUE_PATH}/Raw_reads/ --run-list ${CATALOGUE_PATH}/runs.tsv" | awk '{print $NF}')
  echo "Submitted FetchReadsSlurm $fetch_reads_job_id"
}

function GenerateSamplesheet {
  if [ -n "$fetch_reads_job_id" ]; then
    export dependency_condition="--dependency=afterok:${fetch_reads_job_id}"
  else
    export dependency_condition=""
  fi
  sbatch ${dependency_condition} \
  -J samplesheet_$ASSEMBLY_STUDY_ACC \
  --time=1:00:00 \
  --mem=5G \
  --mail-user=$USER@ebi.ac.uk \
  --mail-type=ALL \
  --wrap="bash ${REPO_PATH}/scripts/run_samplesheet_generation.sh -a $ASSEMBLY_STUDY_ACC -c $CATALOGUE_PATH -r $READS_STUDY_ACC"
}

ASSEMBLY_STUDY_ACC=""
READS_STUDY_ACC=""
CATALOGUE_PATH_INPUT=""
SKIP_FETCH="false" # By default fetch step included
SCRIPT_PATH="$(readlink -f $0)"
REPO_PATH="$(dirname $SCRIPT_PATH)"

while getopts 'a:r:c:f:' flag; do
    case "${flag}" in
        a) ASSEMBLY_STUDY_ACC="$OPTARG" ;;
        r) READS_STUDY_ACC="$OPTARG" ;;
        c) CATALOGUE_PATH_INPUT="$OPTARG" ;;
        f) SKIP_FETCH="$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

CATALOGUE_PATH=$(readlink -f "$CATALOGUE_PATH_INPUT")

if [[ -z $ASSEMBLY_STUDY_ACC || -z $READS_STUDY_ACC || -z $CATALOGUE_PATH || -z $REPO_PATH ]]; then
    echo "Missing required options"
    exit 1
fi


mkdir -p "$CATALOGUE_PATH"
if [[ $SKIP_FETCH == 'false' ]]; then
  . "/hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh"
  mitload fetchtool
  FetchAssembliesSlurm
  GenerateRename
  FetchReadsSlurm
fi

GenerateSamplesheet