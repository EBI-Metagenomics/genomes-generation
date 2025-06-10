process PREPARE_TSV_FOR_UPLOADER {

    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:0.24.1':
        'quay.io/biocontainers/pandas:0.24.1' }"

    input:
    path genomes_euks
    path genomes_proks
    path assembly_software_file
    path stats_euks
    path stats_proks
    path coverage_euks
    path coverage_proks
    path rna
    path taxonomy_euks
    path taxonomy_proks

    output:
    path "*.tsv"                   , optional: true, emit: tsv_for_uploader
    path "unclassified_genomes.txt", optional: true, emit: unclassified_genomes_file
    path "versions.yml"                            , emit: versions

    script:
    def args_genomes_euks = genomes_euks ? "--mags-euks ${genomes_euks}" : "" ;
    def args_genomes_proks = genomes_proks ? "--mags-proks ${genomes_proks}": "" ;
    def args_stats_euks = stats_euks ? "--stats-euks ${stats_euks}": "" ;
    def args_stats_proks = stats_proks ? "--stats-proks ${stats_proks}": "" ;
    def args_coverage_euks = coverage_euks ? "--coverage-euks ${coverage_euks}": "" ;
    def args_coverage_proks = coverage_proks ? "--coverage-proks ${coverage_proks}": "" ;
    def args_assembly_file = assembly_software_file ? "--assembly-software-file ${assembly_software_file}": "" ;
    def args_rna = rna ? "--rna-outs ${rna}": "" ;
    def args_tax_euks = taxonomy_euks ? "--tax-euks ${taxonomy_euks}": "" ;
    def args_tax_proks = taxonomy_proks ? "--tax-proks ${taxonomy_proks}": "" ;
    """
    tsv_for_genome_upload.py \
        $args_genomes_euks \
        $args_genomes_proks \
        $args_stats_euks \
        $args_stats_proks \
        $args_coverage_euks \
        $args_coverage_proks \
        $args_assembly_file \
        ${args_rna} \
        ${args_tax_euks} \
        ${args_tax_proks} \
        --metagenome '$params.metagenome' \
        --biomes '${params.biomes}' \
        --absolute-path "$params.outdir"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
