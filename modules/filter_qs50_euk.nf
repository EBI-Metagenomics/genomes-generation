process FILTER_QS50 {

    publishDir(
        path: "${params.outdir}/qs50",
        mode: "copy"
    )

    label 'process_light'

    input:
    tuple val(name), path(quality_file)
    tuple val(name), path(concoct_bins)
    tuple val(name), path(metabat_bins)
    tuple val(name), path(concoct_bins_merged)
    tuple val(name), path(metabat_bins_merged)

    output:
    tuple val(name), path("output_genomes/*"), emit: qs50_filtered_genomes

    script:
    """
    # prepare drep quality file
    cat ${quality_file} |\
        awk '{{if($2 - 5*$3 >=50){{print $0}}}}' |\
        sort -k 2,3 -n |\
        tr '\\t' ',' > quality_file.csv
    if [[ $(wc -l <quality_file.csv) -lt 2  ]];
    then
        mkdir -p output_genomes
        touch quality_file.csv
    else
        mv ${concoct_bins} ${metabat_bins} ${concoct_bins_merged} ${metabat_bins_merged} output_genomes
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
    fi
    """
}
