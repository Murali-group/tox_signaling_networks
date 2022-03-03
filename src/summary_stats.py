
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from src import toxcast_utils as t_utils
from src import toxcast_settings as t_settings
from src.utils import file_utils as utils


def get_summary_stats(version="2018_01-toxcast-d2d-p1_5-u1_25", 
        summary_file="network_summaries.csv", scope="permute-dir-undir", forced=False):
    """ Function to aggregate summary statistics for every network
    returns a dataframe containing the counted metrics for each chemical
    """
    TOXCAST_DATA = t_utils.loadToxcastData(t_settings.INTERACTOMES[version])
    #inputs_dir = "inputs/%s/" % (version)
    t_settings.set_version(version)
    inputs_dir = t_settings.INPUTSPREFIX 
    outputs_dir = "outputs/%s/weighted" % (version)
    chemicals = utils.readItemList("%s/chemicals.txt" % (inputs_dir), 1)
    #hits_template = "%s/hit-prots/%%s-hit-prots.txt" % (inputs_dir)
    #nonhits_template = "%s/hit-prots/%%s-nonhit-prots.txt" % (inputs_dir)
    #rec_tfs_template = "%s/rec-tfs/%%s-rec-tfs.txt" % (inputs_dir)
    chem_rec, chem_tfs = TOXCAST_DATA.chemical_rec, TOXCAST_DATA.chemical_tfs
    chem_prot_hit_vals = TOXCAST_DATA.chemical_protein_hit
    paths_dir = "%s/edgelinker" % (outputs_dir)
    paths_template = "%s/%%s-paths.txt" % (paths_dir)

    out_dir = "%s/stats/summary" % outputs_dir
    t_utils.checkDir(out_dir)
    summary_file = "%s/%s" % (out_dir, summary_file)
    if os.path.isfile(summary_file) and not forced:
        print("Reading network summary stats from '%s'. Set forced to True to overwrite it." % (summary_file))
        df = pd.read_csv(summary_file, index_col=0) 
    else:
        print("Reading in the stats from the response networks in", paths_dir)
        chemical_names, chemical_name_to_id = t_utils.getChemicalNameMaps()
        chemical_names = {chemical: chemical_names[chemical] for chemical in chemicals}
        chemical_prots = {}
        chemical_num_paths = {}
        chemical_num_edges = {}
        chemical_avg_path_lengths = {}
        chemical_rec = {}
        chemical_tfs = {}
        chemical_net_rec = {}
        chemical_net_tfs = {}
        chemical_hits = {}
        chemical_nonhits = {}
        chemical_net_hits = {}
        chemical_net_nonhits = {}
        chemical_inter_hits = {}
        chemical_inter_nonhits = {}
        chemical_inter_net_hits = {}
        chemical_inter_net_nonhits = {}
        # also get the q-value for each chemical
        chemical_pvals = {} 
        pvals_file = "%s/stats/stat-sig-%s/gpd-pval.txt" % (outputs_dir, scope)
        if os.path.isfile(pvals_file):
            with open(pvals_file, 'r') as file_handle:
                header = file_handle.readline().rstrip().split('\t')
            pval_col = header.index("200") + 1
            chemical_pvals = {chem: pval for chem, pval in utils.readColumns(pvals_file, 1, pval_col)}
        chemical_qvals = {} 
        qvals_file = "%s/stats/stat-sig-%s/bfcorr_pval_qval.txt" % (outputs_dir, scope)
        if os.path.isfile(qvals_file):
            chemical_qvals = t_utils.getPvals(outputs_dir, scope, sig_cutoff_type="FDR")
        for chemical in tqdm(chemicals):
            #prots, paths = getProteins(paths=paths_template % chemical, max_k=200, ties=True)
            paths = t_utils.getPaths(paths_template % chemical, max_k=200, ties=True)
            prots = set()
            num_paths = len(paths)
            edges = set()
            path_lengths = []
            for path in paths:
                path = path.split('|')
                # path length is the number of edges in a path
                path_lengths.append(len(path)-1)
                prots = prots.union(set(path))
                for i in range(len(path)-1):
                    edges.add((path[i], path[i+1]))

            chemical_prots[chemical] = len(prots)
            chemical_num_paths[chemical] = len(paths) 
            chemical_avg_path_lengths[chemical] = np.mean(path_lengths)
            chemical_num_edges[chemical] = len(edges)
            #rec, tfs = t_utils.getRecTFs(rec_tfs_template % chemical)
            rec, tfs = chem_rec[chemical], chem_tfs[chemical]
            chemical_rec[chemical] = len(rec)
            chemical_tfs[chemical] = len(tfs)
            chemical_net_rec[chemical] = len(prots.intersection(rec))
            chemical_net_tfs[chemical] = len(prots.intersection(tfs))
            # read the hits and nonhits for each chemical to calculate how many of them are in the network
            #hits = utils.readItemSet(hits_template % chemical, 1)
            #nonhits = utils.readItemSet(nonhits_template % chemical, 1)
            hits = set([p for p, hit_val in chem_prot_hit_vals[chemical].items() \
                    if hit_val == 1])
            nonhits = set([p for p, hit_val in chem_prot_hit_vals[chemical].items() \
                    if hit_val == 0])
            chemical_hits[chemical] = len(hits)
            chemical_nonhits[chemical] = len(nonhits)
            chemical_net_hits[chemical] = len(hits.intersection(prots))
            chemical_net_nonhits[chemical] = len(nonhits.intersection(prots))
            # subtract the rec and tfs to get just the intermediate hits and nonhits
            chemical_inter_hits[chemical] = len(hits.difference(rec.union(tfs)))
            chemical_inter_nonhits[chemical] = len(nonhits.difference(rec.union(tfs)))
            chemical_inter_net_hits[chemical] = len(hits.intersection(prots).difference(rec.union(tfs)))
            chemical_inter_net_nonhits[chemical] = len(nonhits.intersection(prots).difference(rec.union(tfs)))

        # write these metrics to a file
        df = pd.DataFrame({
            "name": chemical_names,
            "prots": chemical_prots, "num_paths": chemical_num_paths, "pvals": chemical_pvals, "qvals": chemical_qvals,
            "num_edges": chemical_num_edges, "avg_path_lengths": chemical_avg_path_lengths,
            "net_rec": chemical_net_rec, "net_tfs": chemical_net_tfs, "hit_rec": chemical_rec, "hit_tfs": chemical_tfs,
            "net_hits": chemical_net_hits, "net_nonhits": chemical_net_nonhits, 'hits': chemical_hits, 'nonhits': chemical_nonhits,
            "inter_net_hits": chemical_inter_net_hits, "inter_net_nonhits": chemical_inter_net_nonhits, "inter_hits": chemical_inter_hits, "inter_nonhits": chemical_inter_nonhits, 
        })
        print("Writing: ", summary_file)
        df.to_csv(summary_file, header=True, columns=[
            'name', 'prots', 'num_paths', 'num_edges', 'avg_path_lengths', 'hits', 'nonhits', 'net_hits', 'net_nonhits', 'hit_rec', 'hit_tfs', 'net_rec', 'net_tfs', 
            'inter_net_hits', 'inter_net_nonhits', 'inter_hits', 'inter_nonhits', 'pvals', 'qvals'
        ])

    # change the index or chemical id to unicode (string)
    #df.index = df.index.map(unicode) 

    return df

