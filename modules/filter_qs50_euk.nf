process FILTER_QS50 {

    publishDir(
        path: "${params.outdir}/qs50",
        mode: "copy"
    )

    label 'process_light'

    input:
    path concoct_quality
    path metabat_quality
    path concoct_bins
    path metabat_bins

    output:
    path "output_genomes/*", emit: qs50_filtered_genomes


    script:
    """
    # prepare drep quality file
    echo "genome,completeness,contamination" > quality_file.csv
    cat ${concoct_quality}  ${metabat_quality}|\
        awk '{{if($2 - 5*$3 >=50){{print $0}}}}' |\
        sort -k 2,3 -n |\
        tr '\\t' ',' >> quality_file.csv
    if [[ $(wc -l <quality_file.csv) -lt 2  ]];
    then
        touch quality_file.csv
        exit 0
    fi

    mv ${concoct_bins} ${metabat_bins} output_genomes
    ls output_genomes > genomes.txt

    cat quality_file.csv | grep -v "genome,completeness,contaminatio" |\
        cut -f 1 -d , > .f
    grep -f .f genomes.txt | sed -e 's/^/input_genomes\//'  > .g
    mv .g genomes.txt

    ls output_genomes > all_genomes.txt
    grep -v -f genomes.txt output_genomes > remove.txt
    for i in $(cat remove.txt); do
        rm output_genomes/${i};
    done
    """
}
