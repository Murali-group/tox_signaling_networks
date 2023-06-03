# Toxicant Signaling Networks
Code and datasets for building and evaluating toxicant signaling networks.

## Setup
```
conda create -y -n tox_sig_nets python=3.7
conda activate tox_sig_nets
conda install --file requirements.txt
git clone https://github.com/jlaw9/gpdPerm.git
```
To run EdgeLinker, you must have a suitable java installation (> 1.8)

## Supplementary Files:
### EdgeLinker computed paths 
[`outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/edgelinker-outputs.tar.gz`](https://github.com/Murali-group/tox_signaling_networks/blob/main/outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/edgelinker-outputs.tar.gz)
### Interactome
- Interactome: [`inputs/2018_01-toxcast-net/2019-02-18-human-ppi-d2d-dir-weighted-cap0_99.txt.gz`](https://github.com/Murali-group/tox_signaling_networks/blob/main/inputs/2018_01-toxcast-net/2019-02-18-human-ppi-d2d-dir-weighted-cap0_99.txt.gz)
- Iteraction evidence: [`inputs/2018_01-toxcast-net/2018_01interactome-evidence.tsv.gz`](https://github.com/Murali-group/tox_signaling_networks/blob/main/inputs/2018_01-toxcast-net/2018_01interactome-evidence.tsv.gz)
- Confidence scores for each evidence type with > 100 interactions: [`inputs/2018_01-toxcast-net/2019-02-18-edge-source-weights.tsv`](https://github.com/Murali-group/tox_signaling_networks/blob/main/inputs/2018_01-toxcast-net/2019-02-18-edge-source-weights.tsv)
- UniProt to gene mapping: [`inputs/2018_01-toxcast-net/2018_01_uniprot_mapping.tsv`](https://github.com/Murali-group/tox_signaling_networks/blob/main/inputs/2018_01-toxcast-net/2018_01_uniprot_mapping.tsv)
- Chemical identifier (EPA ACToR code) to name mapping: [`inputs/toxcast-tox21-v3/Chemical_Summary_190708.csv`](https://github.com/Murali-group/tox_signaling_networks/blob/main/inputs/toxcast-tox21-v3/Chemical_Summary_190708.csv)
   - e.g., lovastatin: C75330755; BPA: C80057
### Enriched Gene Ontology (GO) terms
- All chemicals: [`outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/stats/go-analysis-k150/389chemicals-sig-terms-bonferroni-c0_01.tsv`](https://github.com/Murali-group/tox_signaling_networks/blob/main/outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/stats/go-analysis-k150/389chemicals-sig-terms-bonferroni-c0_01.tsv)

## Regenerate Results
The commands below outline the steps necessary to regenerate the results and figures from the paper.

### Create Toxicant Signaling Networks with EdgeLinker
```
python masterscript.py --version 2018_01-toxcast-d2d-p1_5-u1_25 --response-networks
```

### Network Summary Plot
- To generate Figure 1a-d:
```
python src/plot/version_plots/network_summary_plots.py  --version 2018_01-toxcast-d2d-p1_5-u1_25 -c --pdf  -k 150
```

### Statistical Significance - Randomized interactome
- To generate the "random networks", the command below creates 1,000 randomized interactomes using edge swapping
- Then run EdgeLinker, using the original set of sources and targets for each toxicant, on each of those randomized interactomes
- The `--super-computer` option creates many PBS scripts and submits them to a cluster at VT (baobab)
```
version="2018_01-toxcast-d2d-p1_5-u1_25"
python masterscript.py \
    --version $version \
    --random-networks \
    --scope permute-dir-undir \
    --num-random-networks 1 1000  \
    --parallel-processes 10   \
    --super-computer \
    --edgelinker-k 1000 \
    --cleanup \
    -k 50 -k 100 -k 150 -k 200 \
    --forcealg
```
- This is the command to compute the p-values after the random networks have been generated
```
python masterscript.py --version $version --stat-sig -k 50 -k 100 -k 150 -k 200 --scope permute-dir-undir
```
- outputs will be here: `outputs/2018_01-toxcast-d2d-p3-u1_25/weighted/stats/stat-sig-permute-dir-undir/bfcorr_pval_qval.txt`
- Use this command to get the chemicals that have a q-value < 0.01
```
python masterscript.py --version $version --scope permute-dir-undir --sig-cutoff 0.01 --sig-cutoff-type FDR
# alternatively:
python src/stat_sig/compute_stat_sig.py  --chemicals inputs/versions/$version/chemicals.txt  --compute-score-pval  --num-random 1 1000   -k 50 -k 100 -k 150 -k 200   --out-dir outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted//stats/stat-sig-permute-dir-undir/ --pval-cutoff 0.01
```

### Overlap Graph
- To generate the graph with nodes and edges present in more than 1/3 of the toxicant signaling networks, use these commands
```
python src/plot/overlap/plot_overlap_heatmap.py  --chemicals inputs/versions/$version//sig-chemicals.txt  --out-dir  outputs/$version/weighted//stats/overlap/permute-dir-undir-FDR-0_05-k150/ --k-limit 150 --paths-dir outputs/$version/weighted//edgelinker/  --write-overlap  --plot-overlap  --plot-hits-overlap --mapping-file inputs/2018_01-toxcast-net/2018_01_uniprot_mapping.tsv

python src/graphspace/post_to_graphspace_overlap.py     --ppi inputs/versions/$version/$version-interactome-undir-penalty.txt     --version $version --ev-file /home/jeffl/svnrepo/data/interactions/compiled/2018_01/with-d2d/2018_01pathlinker-nokegg.tsv  --sourcetarget inputs/versions/$version/rec-tfs/all-rec-tfs.txt     --casestudy      --edge-overlap outputs/$version/weighted/stats/overlap/permute-dir-undir-FDR-0_05-k150/overlapping-edges.txt     --node-overlap outputs/$version/weighted/stats/overlap/permute-dir-undir-FDR-0_05-k150/overlapping-proteins.txt    --ctd-support inputs/ctd/CTD_chem_gene_ixns_human_phospho.tsv   --tag $version  --sig-chemicals inputs/versions/$version//sig-chemicals.txt --mapping-file inputs/2018_01-toxcast-net/2018_01_uniprot_mapping.tsv --overlap-cutoff 0.33 --graph-name overlap-graph-sig-0_33-k150
```

### CTD Overlap
```
version="2018_01-toxcast-d2d-p1_5-u1_25"; python src/ctd_analysis/run_ctd_analysis.py \
    --chemicals inputs/versions/$version/sig-chemicals.txt \
    --out-dir outputs/$version/weighted/stats/ctd-analysis-test/ \
    --paths outputs/$version/weighted/edgelinker/ \
    --k-limit 150 \
    --ctd-file inputs/ctd/CTD_chem_gene_ixns_human_phospho.tsv \
    --interactome inputs/2018_01-toxcast-net/2019-02-18-human-ppi-weighted-cap0_99.txt  \
    --mapping-file inputs/2018_01-toxcast-net/2018_01_uniprot_mapping.tsv
```

### GO term enrichment
- all chemicals:
```
version="2018_01-toxcast-d2d-p1_5-u1_25"; 
python src/goterm_analysis/run_david_analysis.py \
  --chemicals inputs/versions/$version/chemicals.txt \
  --out-dir outputs/$version/weighted/stats/go-analysis-k150/ \
  --paths outputs/$version/weighted/edgelinker/ \
  --k-limit 150 \
  --correction-type BF
```
- For significant chemicals only: change the `--chemicals` file to `sig-chemicals.txt`
- Use `--pval-cutoff 0.01` to set the BF-corrected p-value cutoff to 0.01

### Case Studies
#### Settings used to run REVIGO
- We set `allowed semantic similarity` to 0.5, which corresponds to a "small" list size
- For the database of GO term sizes, we used the `homo sapiens` database
- All other settings set to defaults
##### choosing the terms after revigo
- We removed terms with a frequency > 5%
- Many terms still represented similar functions, so we manually chose a single term to represent clusters in the revigo output.
#### GO term heatmaps
- See the notebook `src/jupyter-notebooks/goterm-heatmap.ipynb`

#### Post to GraphSpace
##### Lovastatin
- The ID for lovastatin is `C75330755`
```
python src/graphspace/post_to_graphspace_wrapper.py  \
    --revigo-file outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/stats/go-analysis-k150/revigo/lovastatin/lovastatin-selected-terms.tsv \
    --ctd-support-file=inputs/ctd/CTD_chem_gene_ixns_human_phospho.tsv  \
    --version 2018_01-toxcast-d2d-p1_5-u1_25 -S C75330755 \
    --user <username> \
    --pass <password> \
    --parent-nodes \
    --case-study  \
    --postfix=-revigo-ctd
```
##### BPA
- The ID for BPA is `C80057`
```
python src/graphspace/post_to_graphspace_wrapper.py  \
    --revigo-file outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/stats/go-analysis-k150/revigo/bpa/bpa-selected-terms.tsv \
    --ctd-support-file=inputs/ctd/CTD_chem_gene_ixns_human_phospho.tsv  \
    --version 2018_01-toxcast-d2d-p1_5-u1_25 -S C80057 \
    --user <username> \
    --pass <password> \
    --parent-nodes \
    --case-study  \
    --postfix=-revigo-ctd
```
### Post Networks to GraphSpace
- First, make an account on graphspace as well as a group (optional)
- The bash script below is an example of how to call `post_to_graphspace_wrapper.py` for each chemical
```
bash src/graphspace/post-to-gs.sh
```
