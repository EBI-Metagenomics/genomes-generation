#!/usr/bin/env bash

function FetchAssembliesAndReads {
  echo 'Fetching reads and assemblies...'
  . "/hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh"
  mitload fetchtool
  bsub -I -q production "bash fetch-assemblies-tool.sh -v -p $SAMPLE -d ${CATALOGUE_PATH}/Assemblies/"
  bsub -I -q production "bash fetch-reads-tool.sh -v -p $READS_ACC -d ${CATALOGUE_PATH}/Raw_reads/"
}


function Unzip {
  echo 'Unzipping...'
  gunzip -f ${CATALOGUE_PATH}/Assemblies/${SAMPLE}/raw/*gz
}


function RunRenamingScript {
  echo 'Starting renaming script...'
  . /hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh
  mitload assembly_pipeline
  sleep 1
  mkdir -p ${CATALOGUE_PATH}/Uploaded_Assembly_IDs

  python3 ${REPO_PATH}/download_data/rename-erz.py \
  -d ${CATALOGUE_PATH}/Assemblies/${SAMPLE}/raw/ -o ${CATALOGUE_PATH}/Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt

  cat ${CATALOGUE_PATH}/Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt | tr ',' '\t' > ${CATALOGUE_PATH}/rename.tsv
  export CONVERT=${CATALOGUE_PATH}/Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt
}


function Rename {
  while read line; do
    if [[ $line == [SED]RR* ]]
    then
      export OLD=$(echo $line | cut -d ',' -f1)
      export NEW=$(echo $line | cut -d ',' -f2)
      mv ${CATALOGUE_PATH}/Assemblies/${SAMPLE}/raw/${NEW}.fasta ${CATALOGUE_PATH}/Assemblies/${SAMPLE}/raw/${OLD}.fasta
    fi
  done < $CONVERT
}


#function ChangeERRtoERZinReads {
#  for read_file in (ls $CATALOGUE/Raw_reads/${SAMPLE}/raw/* )
#    grep "${name}" ${rename_file} > help_file
#    export from_accession=\$(cat help_file | cut -f1)
#    export to_accession=\$(cat help_file | cut -f2)
#    echo "\${from_accession} --> \${to_accession}"
#
#    zcat "${input_ch[0]}" | sed "s/\${from_accession}/\${to_accession}/g" | gzip > ${run_accession}_changed.fastq.gz
#}
SAMPLE=""
READS_ACC=""
CATALOGUE_PATH=""
SKIP_FETCH="true" # Default to true
REPO_PATH=""

while getopts 'a:r:c:fp:' flag; do
    case "${flag}" in
        a) SAMPLE="$OPTARG" ;;
        r) READS_ACC="$OPTARG" ;;
        c) CATALOGUE_PATH="$OPTARG" ;;
        f) SKIP_FETCH='false' ;;
        p) REPO_PATH="$OPTARG" ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

if [[ -z $SAMPLE || -z $READS_ACC || -z $CATALOGUE_PATH || -z $REPO_PATH ]]; then
    echo "Missing required options"
    exit 1
fi

mkdir -p "$CATALOGUE_PATH"

if [[ $SKIP_FETCH == 'false' ]]; then
  FetchAssembliesAndReads
  Unzip
fi

RunRenamingScript
Rename
#ChangeERRtoERZinReads