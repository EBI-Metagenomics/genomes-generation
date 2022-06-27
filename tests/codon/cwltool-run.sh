#!/bin/bash

set -e

. /hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/mitrc.sh

mitload miniconda

module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
conda activate toil-5.6.0

export OUTDIR=result

while getopts :o:n:y:c:p:d: option; do
	case "${option}" in
		o) OUTDIR=${OPTARG};;
		n) NAME=${OPTARG};;
		y) YML=${OPTARG};;
		c) CWL=${OPTARG};;
		p) MAIN_PATH=${OPTARG};;
	    d) DEBUG=${OPTARG};;
	esac
done

if [ "${DEBUG}" == "True" ]; then
    cwltool --singularity --preserve-entire-environment --debug --leave-container --outdir ${OUTDIR}/${NAME} ${CWL} ${YML}
else
    cwltool --singularity --preserve-entire-environment --leave-container --outdir ${OUTDIR}/${NAME} ${CWL} ${YML}
fi