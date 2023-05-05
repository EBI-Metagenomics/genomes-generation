checkm
parse_checkm

echo "genome,completeness,contamination" > all.stats.clean
cut -f1-3 ../checkm_results.tab | sed 's/\s/,/g' | sed 's/\,/\.fa\,/' >> all.stats.clean