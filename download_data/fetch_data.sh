#!/usr/bin/env bash

function FetchAssembliesAndReads {
  echo 'Fetching reads and assemblies...'
  . "/hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh"
  mitload fetchtool
  bsub -I -q production "bash fetch-assemblies-tool.sh -v -p $SAMPLE -d Assemblies/"
  bsub -I -q production "bash fetch-reads-tool.sh -v -p $READS_ACC -d Raw_reads/"
}


function Unzip {
  echo 'Unzipping...'
  gunzip -f Assemblies/${SAMPLE}/raw/*gz
}


function RunRenamingScript {
  echo 'Starting renaming script...'
  mitload assembly_pipeline

  python3 rename-erz.py \
  -d Assemblies/${SAMPLE}/raw/ -o Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt

  export CONVERT=Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt

  if [[ ! -f $CONVERT ]]; then
    echo 'ERZ to run accession conversion file was not generated successfully'
    exit 1
  fi
}


function Rename {
  while read line; do
    if [[ $line == [SED]RR* ]]
    then
      export OLD=$(echo $line | cut -d ',' -f1)
      export NEW=$(echo $line | cut -d ',' -f2)
      mv Assemblies/${SAMPLE}/raw/${NEW}.fasta Assemblies/${SAMPLE}/raw/${OLD}.fasta
    fi
  done < $CONVERT
}


while getopts 'a:r:' flag; do
    case "${flag}" in
        a) export SAMPLE=$OPTARG ;;
        r) export READS_ACC=$OPTARG ;;
    esac
done

FetchAssembliesAndReads
Unzip
RunRenamingScript
Rename