# Toxicant Signaling Networks
Code and datasets for building and evaluating toxicant signaling networks.

## Setup
```
conda create -y -n tox_sig_nets python=3.7
conda activate tox_sig_nets
conda install --file requirements.txt
```
To run EdgeLinker, you must have a suitable java installation (> 1.8)

## Regenerate Results
The commands below outline the steps necessary to regenerate the results and figures from the paper.

### Create Toxicant Signaling Networks with EdgeLinker
```
python masterscript.py --version 2018_01-toxcast-d2d-p1_5-u1_25 --response-networks
```
### Statistical Significance - Randomized interactome
     - To generate the "random networks", the command below creates 1,000 randomized interactomes using edge swapping
     - Then run EdgeLinker, using the original set of sources and targets for each toxicant, on each of those randomized interactomes
     - The `--super-computer` option creates many PBS scripts and submits them to a cluster at VT (baobab)
```
version="2018_01-toxcast-d2d-p1_5-u1_25"
python master-script.py \
    --version $version \
    --random-networks \
    --scope permute-dir-undir \
    --num-random-networks 1 1000  \
    --parallel-processes 10   \
    --super-computer \
    --edgelinker-k 1000 \
    --cleanup \
    -k 50 -k 100 -k 200 \
    --forcealg
```
     - This is the command to compute the p-values after the random networks have been generated
```
python master-script.py --version $version --stat-sig -k 50 -k 100 -k 200 --scope permute-dir-undir
```
     - outputs will be here: `outputs/2018_01-toxcast-d2d-p3-u1_25/weighted/stats/stat-sig-permute-dir-undir/bfcorr_pval_qval.txt`
### CTD Overlap
```
version="2018_01-toxcast-d2d-p1_5-u1_25"; python src/ctd_analysis/run_ctd_analysis.py \
    --chemicals inputs/versions/$version/sig-chemicals.txt \
    --out-dir outputs/$version/weighted/stats/ctd-analysis-test/ \
    --paths outputs/$version/weighted/edgelinker/ \
    --ctd-file inputs/ctd/2020-03-02/CTD_chem_gene_ixns.tsv \
    --interactome inputs/2018_01-toxcast-net/2019-02-18-human-ppi-weighted-cap0_99.txt  \
    --mapping-file inputs/2018_01-toxcast-net/2018_01_uniprot_mapping.tsv
```
### GO term enrichment
     - all chemicals:
```
version="2018_01-toxcast-d2d-p1_5-u1_25"; 
python src/goterm_analysis/run_david_analysis.py \
  --chemicals inputs/versions/$version/chemicals.txt \
  --out-dir outputs/$version/weighted/stats/go-analysis/ \
  --paths outputs/$version/weighted/edgelinker/ \
  --correction-type BF
```
     - For significant chemicals only: add the option `--pval-cutoff 0.01`
### Fig 2e - Overlap Graph
#### TODO Post to GraphSpace
### Case Studies
#### Settings used to run REVIGO
     - We set `allowed semantic similarity` to 0.4, which corresponds to a "tiny" list size
     - For the database of GO term sizes, we used the `homo sapiens` database
     - All other settings set to defaults
##### choosing the terms after revigo
      - We removed terms with a frequency > 5%
      - Many terms still represented similar functions, so we manually chose a single term to represent clusters in the revigo output.
#### Post to GraphSpace
##### Lovastatin
      - The ID for lovastatin is `C75330755`
```
python src/graphspace/post_to_graphspace_wrapper.py  \
    --revigo-file outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/stats/go-analysis/revigo/lovastatin/revigo-lovastatin-selected-terms.csv \
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
    --revigo-file outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/stats/go-analysis/revigo/bpa/bpa-revigo-selected-terms.csv \
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
```
python src/graphspace/post_to_graphspace_wrapper.py \
    --ctd-support-file=inputs/ctd/CTD_chem_gene_ixns_human_phospho.tsv \
    --version 2018_01-toxcast-d2d-p1_5-u1_25 \
    -S C137304 \
    --user <username> --pass <password> \
    --parent-nodes --case-study \
    --term-counts-file=outputs/2018_01-toxcast-d2d-p1_5-u1_25/weighted/stats/go-analysis/505chemicals-sig-terms-bonferroni-c0_01-counts.tsv \
    --group=<group-name>
```
