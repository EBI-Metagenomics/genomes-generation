/*
 * Predict bacterial 5S, 16S, 23S rRNA and tRNA genes
*/
process DETECT_RRNA {

    tag "${fasta}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.detect_rrna:v3.1'

    input:
    path(fasta)
    path cm_models

    output:
    path('results_folder/*'), emit: rrna_out_results

    script:
    """
    shopt -s extglob

    RESULTS_FOLDER=results_folder
    FASTA=${fasta}
    CM_DB=${cm_models}

    BASENAME=\$(basename "\${FASTA}")
    FILENAME="\${BASENAME%.*}"

    mkdir -p "\${RESULTS_FOLDER}"

    echo "[ Detecting rRNAs ] "

    for CM_FILE in "\${CM_DB}"/*.cm; do
        MODEL=\$(basename "\${CM_FILE}")
        echo "Running cmsearch for \${MODEL}..."
        cmsearch -Z 1000 \
            --hmmonly \
            --cut_ga --cpu ${task.cpus} \
            --noali \
            --tblout "\${RESULTS_FOLDER}/\${FILENAME}_\${MODEL}.tblout" \
            "\${CM_FILE}" "\${FASTA}" 1> "\${RESULTS_FOLDER}/\${FILENAME}_\${MODEL}.out"
    done

    echo "Concatenating results..."
    cat "\${RESULTS_FOLDER}/\${FILENAME}"_*.tblout > "\${RESULTS_FOLDER}/\${FILENAME}.tblout"

    echo "Removing overlaps..."
    cmsearch-deoverlap.pl \
    --maxkeep \
    --clanin "\${CM_DB}/ribo.claninfo" \
    "\${RESULTS_FOLDER}/\${FILENAME}.tblout"

    mv "\${FILENAME}.tblout.deoverlapped" "\${RESULTS_FOLDER}/\${FILENAME}.tblout.deoverlapped"

    echo "Parsing final results..."
    parse_rRNA-bacteria.py -i \
    "\${RESULTS_FOLDER}/\${FILENAME}.tblout.deoverlapped" 1> "\${RESULTS_FOLDER}/\${FILENAME}_rRNAs.out"

    rRNA2seq.py -d \
    "\${RESULTS_FOLDER}/\${FILENAME}.tblout.deoverlapped" \
    -i "\${FASTA}" 1> "\${RESULTS_FOLDER}/\${FILENAME}_rRNAs.fasta"

    echo "[ Detecting tRNAs ]"
    tRNAscan-SE -B -Q \
    -m "\${RESULTS_FOLDER}/\${FILENAME}_stats.out" \
    -o "\${RESULTS_FOLDER}/\${FILENAME}_trna.out" "\${FASTA}"

    parse_tRNA.py -i "\${RESULTS_FOLDER}/\${FILENAME}_stats.out" 1> "\${RESULTS_FOLDER}/\${FILENAME}_tRNA_20aa.out"

    echo "Completed"
    """

    stub:
    """
    mkdir results_folder
    cd results_folder
    touch ACC_bin.1_clean_RF00001.cm.out \
        ACC_bin.1_clean_RF00177.cm.out \
        ACC_bin.1_clean_RF02541.cm.out \
        ACC_bin.1_clean_rRNAs.fasta \
        ACC_bin.1_clean_rRNAs.out \
        ACC_bin.1_clean_stats.out \
        ACC_bin.1_clean.tblout.deoverlapped \
        ACC_bin.1_clean_tRNA_20aa.out \
        ACC_bin.1_clean_trna.out
    """
}
