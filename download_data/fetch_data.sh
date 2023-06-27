#!/usr/bin/env bash

function FetchAssembliesAndReads {
  echo 'Fetching reads and assemblies...'
  . "/hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh"
  mitload fetchtool
  bsub -I -q production "bash fetch-assemblies-tool.sh -v -p $SAMPLE -d $CATALOGUE/Assemblies/"
  bsub -I -q production "bash fetch-reads-tool.sh -v -p $READS_ACC -d $CATALOGUE/Raw_reads/"
}


function Unzip {
  echo 'Unzipping...'
  gunzip -f $CATALOGUE/Assemblies/${SAMPLE}/raw/*gz
}


function RunRenamingScript {
  echo 'Starting renaming script...'
  . /hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh
  mitload assembly_pipeline
  mkdir -p $CATALOGUE/Uploaded_Assembly_IDs

  python3 ./download_data/rename-erz.py \
  -d $CATALOGUE/Assemblies/${SAMPLE}/raw/ -o $CATALOGUE/Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt

  cat $CATALOGUE/Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt | tr ',' '\t' > $CATALOGUE/rename.tsv
  export CONVERT=$CATALOGUE/Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt

  if [[ ! -f $CONVERT ]]; then
    echo 'ERZ to run accession conversion file was not generated successfully'
    exit 1
  fi
  rm -rf $CATALOGUE/Uploaded_Assembly_IDs
}


function Rename {
  while read line; do
    if [[ $line == [SED]RR* ]]
    then
      export OLD=$(echo $line | cut -d ',' -f1)
      export NEW=$(echo $line | cut -d ',' -f2)
      mv $CATALOGUE/Assemblies/${SAMPLE}/raw/${NEW}.fasta $CATALOGUE/Assemblies/${SAMPLE}/raw/${OLD}.fasta
    fi
  done < $CONVERT
}


while getopts 'a:r:c:f:' flag; do
    case "${flag}" in
        a) export SAMPLE=$OPTARG ;;
        r) export READS_ACC=$OPTARG ;;
        c) export CATALOGUE=$OPTARG ;;
        f) SKIP_FETCH='true' ;;
    esac
done

mkdir -p $CATALOGUE
if [[ $SKIP_FETCH = false ]]
then
  FetchAssembliesAndReads
  Unzip
fi
RunRenamingScript
Rename