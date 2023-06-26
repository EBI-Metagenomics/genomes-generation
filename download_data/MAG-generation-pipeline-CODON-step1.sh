#!/usr/bin/env bash

# The script downloads all raw data for bin and MAG generation, renames files and cleans reads
# Run as emgpr with logname emgpr
# the directory where the script is executed needs to contain the following subdirectories:
# Assemblies
# Raw_reads
# Logs
# Uploaded_Assembly_IDs

# The Uploaded_Assembly_IDs directory will contain uploaded_runs.txt file renamed as
# ${secondary_assembly_study_acc}.uploaded_runs.txt


function Usage {
    echo "Usage: $0 [-a [E|R|D|]XXXXXX] [-r [E|R|D|]XXXXXX] [-h /path/to/host/genome] [-s] [-c] [-f] [-u]"
    echo "Options:"
    echo "-a   Secondary accession of the assembly study"
    echo "-r   Secondary accession of the raw reads study"
    echo "-h   Path to the host genome fasta; must have bwa index"
    echo "-s   Use this flag if the reads are single-end"
    echo "-c   Use this flag if you already have clean reads"
    echo "-f   Use this flag to skip assembly and read fetching"
    echo "-u   Use this flag to skip ERZ to [E|D|S]RR file generation"
    exit 1
}


function GenerateDirectories {
  if [[ ! -d "Logs/${SAMPLE}" ]]
  then
    mkdir Logs/${SAMPLE}
    mkdir Logs/${SAMPLE}/Cleaning
    mkdir Logs/${SAMPLE}/BinCleaning
    mkdir Logs/${SAMPLE}/MetaWRAP
    mkdir Logs/${SAMPLE}/Taxonomy
    mkdir Logs/${SAMPLE}/drep
  fi
}


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
  if [[ $SKIP_CONVERSION_GENERATION = false ]]
  then
    python /hps/software/users/rdf/metagenomics/service-team/repos/mags-generation-pipeline/scripts/rename-erz.py \
    -d Assemblies/${SAMPLE}/raw/ -o Uploaded_Assembly_IDs/${SAMPLE}.uploaded_runs.txt
  fi
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


function CleanFastq {
  echo "Cleaning Fastq..."
  for ACC in $ALL_SAMPLES
  do
    if [[ $SINGLE = true ]] || [[ ! -f Raw_reads/${READS_ACC}/raw/${ACC}'_1.fastq.gz' ]]
    then
      bsub -J "FASTQ-CLEANING-${ACC}" -M 30000 -n 8 -q production -o Logs/${SAMPLE}/Cleaning/${ACC}.log \
      "bash /hps/software/users/rdf/metagenomics/service-team/repos/mags-generation-pipeline/scripts/clean_fastq.sh \
      -t 8 -f Raw_reads/${READS_ACC}/raw/${ACC}'.fastq.gz' -c ${HOST_PATH}"
    else
      bsub -J "FASTQ-CLEANING-${ACC}" -M 30000 -n 8 -q production -o Logs/${SAMPLE}/Cleaning/${ACC}.log \
      "bash /hps/software/users/rdf/metagenomics/service-team/repos/mags-generation-pipeline/scripts/clean_fastq.sh \
      -t 8 -f Raw_reads/${READS_ACC}/raw/${ACC}'_1.fastq.gz' -r Raw_reads/${READS_ACC}/raw/${ACC}'_2.fastq.gz' \
      -c ${HOST_PATH}"
    fi
  done
}


function WaitCleaning {
  echo "Waiting for cleaning jobs to complete..."
  echo "(You can exit out of this script and monitor the completion of the jobs yourself before starting step 2)"
  for ACC in $ALL_SAMPLES
  do
    bwait -w "ended(FASTQ-CLEANING-${ACC})"
    echo "Done waiting on job FASTQ-CLEANING-${ACC}"
  done
}


while getopts 'a:r:h:scfu' flag; do
    case "${flag}" in
        a) export SAMPLE=$OPTARG ;;
        r) export READS_ACC=$OPTARG ;;
        h) export HOST_PATH=$OPTARG ;;
        s) SINGLE='true' ;;
        c) SKIP_CLEANING='true' ;;
        f) SKIP_FETCH='true' ;;
        u) SKIP_CONVERSION_GENERATION='true' ;;
    esac
done

if [[ -z $SAMPLE ]] || [[ -z $READS_ACC ]]; then
  echo 'Not all of the arguments are provided'
  Usage
fi


GenerateDirectories
if [[ $SKIP_FETCH = false ]]
then
  FetchAssembliesAndReads
  Unzip
fi
RunRenamingScript
Rename
if [[ $SKIP_CLEANING = false ]]
then
  export ALL_SAMPLES=$(ls Assemblies/${SAMPLE}/raw/*fasta | cut -d '/' -f4 | cut -d '.' -f1)
  CleanFastq
  WaitCleaning
fi



