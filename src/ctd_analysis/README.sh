#Example call: 
echo """version="2018_01-toxcast-d2d-p1_5-u1_25";
python src/ctd_analysis/run_ctd_analysis.py \
    --chemicals inputs/versions/$version/sig-chemicals.txt \
    --out-dir outputs/$version/weighted/stats/ctd-analysis/ \
    --paths outputs/$version/weighted/edgelinker/ \
    --ctd-file inputs/ctd/2020-03-02/CTD_chem_gene_ixns.tsv.gz \
    --interactome inputs/2018_01-toxcast-net/2019-02-18-human-ppi-weighted-cap0_99.txt  \
    --mapping-file inputs/2018_01-toxcast-net/2018_01_uniprot_mapping.tsv
"""
version="2018_01-toxcast-d2d-p1_5-u1_25";
python src/ctd_analysis/run_ctd_analysis.py \
    --chemicals inputs/versions/$version/sig-chemicals.txt \
    --out-dir outputs/$version/weighted/stats/ctd-analysis/ \
    --paths outputs/$version/weighted/edgelinker/ \
    --ctd-file inputs/ctd/2020-03-02/CTD_chem_gene_ixns.tsv.gz \
    --interactome inputs/2018_01-toxcast-net/2019-02-18-human-ppi-weighted-cap0_99.txt  \
    --mapping-file inputs/2018_01-toxcast-net/2018_01_uniprot_mapping.tsv
