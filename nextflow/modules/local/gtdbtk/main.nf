process GTDBTK {

    container 'quay.io/microbiome-informatics/gtdb-tk:2.1.0'
    containerOptions "--bind ${gtdbtk_refdata}:/opt/gtdbtk_refdata"

    label 'process_bigmem'

    input:
    path(genomes_fna, stageAs: "genomes_dir/*")
    path gtdbtk_refdata

    output:
    path("Taxonomy")

    script:
    // TODO: tweak the cpus based on the number of genomes
    """
    GTDBTK_DATA_PATH=/opt/gtdbtk_refdata \
    gtdbtk classify_wf \
    --cpus ${task.cpus} \
    --pplacer_cpus ${task.cpus} \
    --genome_dir genomes_dir \
    --extension fa \
    --out_dir Taxonomy
    """

    stub:
    """
    mkdir Taxonomy

    mkdir -p Taxonomy/classify
    touch Taxonomy/classify/gtdbtk.bac120.summary.tsv
    touch Taxonomy/classify/gtdbtk.ar122.summary.tsv

    echo "user_genome	classification	fastani_reference	fastani_reference_radius	fastani_taxonomy	fastani_ani	fastani_af	closest_placement_reference	closest_placement_radius	closest_placement_taxonomy	closest_placement_ani	closest_placement_af	pplacer_taxonomy	classification_method	note	other_related_references(genome_id,species_name,radius,ANI,AF)	msa_percent	translation_table	red_value	warnings" > gtdbtk_results/classify/gtdbtk.bac120.summary.tsv
    """
}
