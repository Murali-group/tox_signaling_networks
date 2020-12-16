# 2020-03-02 Jeff L: 
# I downloaded the CTD chemical interactions from here:
wget http://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz

# then I preprocessed the file to get just the Human genes, and interactions
# involving phosphorylation:
grep -P "\t9606\t" CTD_chem_gene_ixns.tsv  > processed/CTD_chem_gene_ixns_human.tsv
grep -i "phosphorylation" processed/CTD_chem_gene_ixns_human.tsv | grep -P  "\tprotein\t" > processed/CTD_chem_gene_ixns_human_phospho.tsv


# 2020-09-23 Jeff L: 
# I downloaded the CTD GO term data from here:
wget http://ctdbase.org/reports/CTD_chem_go_enriched.tsv.gz

