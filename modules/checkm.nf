process CHECKM {

    container 'quay.io/microbiome-informatics/genomes-pipeline.checkm:v1'

    input:
    path bins
    path marker_file
    path checkm_db

    output:
    path "bins_taxonomy.tab", emit: bins_taxonomy
    path "bins_qa.tab", emit: bins_qa

    script:
    """
    echo "checkm tree"
    checkm tree -t 8 ${bins} checkm_output

    echo "checkm tree_qa"
    checkm tree_qa checkm_output --tab_table -f bins_taxonomy.tab

    echo "checkm lineage_set"
    checkm lineage_set checkm_output ${marker_file}

    echo "checkm analyze"
    checkm analyze -t 8 ${marker_file} ${bins} checkm_output

    echo "checkm qa"
    checkm qa ${marker_file} checkm_output --tab_table -f bins_qa.tab
    """
}