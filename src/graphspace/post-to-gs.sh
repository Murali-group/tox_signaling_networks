
# post all of the chemicals
for chem in `cut -f 1 inputs/versions/2018_01-toxcast-d2d-p1_5-u1_25/chemicals.txt` ; do 
    echo "$chem"
    echo """python src/graphspace/post_to_graphspace_wrapper.py  \
    --ctd-support-file=inputs/ctd/CTD_chem_gene_ixns_human_phospho.tsv \
    --version 2018_01-toxcast-d2d-p1_5-u1_25 -S $chem \
    --user <username> --pass <password> \
    --parent-nodes --case-study  \
    --group 2022-toxicant-signaling-networks  \
    --term-counts-file=outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/stats/go-analysis/505chemicals-sig-terms-bonferroni-c0_01-counts.tsv
"""
    python src/graphspace/post_to_graphspace_wrapper.py  \
        --ctd-support-file=inputs/ctd/CTD_chem_gene_ixns_human_phospho.tsv \
        --version 2018_01-toxcast-d2d-p1_5-u1_25 -S $chem \
        --user <username> --pass <password> \
        --parent-nodes --case-study  \
        --group 2022-toxicant-signaling-networks  \
        --term-counts-file=outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/stats/go-analysis/505chemicals-sig-terms-bonferroni-c0_01-counts.tsv
done
echo "DONE"
