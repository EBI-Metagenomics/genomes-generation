#!/bin/bash

usage()
{
cat << EOF
usage: $0 options
Script chooses good quality genomes according to eukcc.csv from concoct and metabat2.
OPTIONS:
   -c      Eukcc concoct.csv [REQUIRED]
   -m      Eukcc metabat2.csv [REQUIRED]
   -a      Directory with concoct genomes [REQUIRED]
   -b      Directory with metabat2 genomes [REQUIRED]
EOF
}

# variables
concoct_csv=
metabat2_csv=
concoct_genomes=
metabat2_genomes=

while getopts “c:m:a:b:” OPTION
do
     case ${OPTION} in
         c)
             concoct_csv=${OPTARG}
             ;;
         m)
             metabat2_csv=${OPTARG}
             ;;
         a)
             concoct_genomes=${OPTARG}
             ;;
         b)
             metabat2_genomes=${OPTARG}
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

out_dir="genomes"
output_quality="quality.csv"
output_genomes="genomes.txt"

mkdir -p ${out_dir} && rm -r  ${out_dir} && mkdir -p ${out_dir}

# prepare drep quality file
echo "genome,completeness,contamination" > ${output_quality}
cat ${concoct_csv}  ${metabat2_csv}|\
    awk '{{if($2 - 5*$3 >=50){{print $0}}}}' |\
    sort -k 2,3 -n |\
    tr '\t' ',' >> ${output_quality}

if [[ $(wc -l <${output_quality}) -lt 2  ]];
then
    touch ${output_quality}
    touch ${output_genomes}
    mkdir -p ${out_dir}/dereplicated_genomes || true
    mkdir -p ${out_dir}/input_bins || true
    exit 0
fi

cat ${output_quality}
mkdir ${out_dir}/input_bins
#ln -s {input.concoct}/*.fa . || true
ln -s ${concoct_genomes}/*.fa ${out_dir}/input_bins/. || true
#ln -s {input.metabat2}/*.fa . || true
ln -s ${metabat2_genomes}/*.fa ${out_dir}/input_bins/. || true

ls ${out_dir}/input_bins/ > genomes.txt

cat ${output_quality} | grep -v "genome,completeness,contamination" | cut -f 1 -d , > .f

grep -f .f genomes.txt > .g
mv .g genomes.txt
rm .f
rm -rf ${out_dir}/input_bins

grep -w -f genomes.txt ${concoct_csv} > concoct_temp.csv
for i in $(cat concoct_temp.csv | cut -f1); do
    cp ${concoct_genomes}/${i} ${out_dir};
done

grep -w -f genomes.txt ${metabat2_csv} > metabat2_temp.csv
for i in $(cat metabat2_temp.csv | cut -f1); do
    cp ${metabat2_genomes}/${i} ${out_dir};
done

rm concoct_temp.csv metabat2_temp.csv