#! /usr/bin/env python

# script to make a histogram of the edge weights and -log of the edge costs
# useful to compare how chaning weights changes the distribution of costs in the -log space

print("Importing libraries")

import os
import sys
import shutil
from optparse import OptionParser
from src.utils import file_utils as utils
from src import toxcast_utils as t_utils
from src import toxcast_settings as t_settings
from math import log
from tqdm import tqdm
# import plotting tools
import matplotlib
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
plt.style.use('ggplot')
#import seaborn as sns
#sns.set_style("whitegrid")
#import numpy as np
# use this to fix the large amount of DEBUG warning messages 
import logging
logging.getLogger('matplotlib.font_manager').disabled = True


parser = OptionParser()
parser.add_option('','--version',type='string',action="append",
        help="Version of the PPI to use. Multiple can be run. Options are: %s." % (', '.join(t_settings.ALLOWEDVERSIONS)))
parser.add_option('','--ev-file',type='string',action="append",
                  help="path/to/evidence_file.txt that contains the edge direction info. If not specified, default is used.")
parser.add_option('-c', '--compare-versions', action='store_true',
                  help='plots will be written to viz/version_plots/edge-weights/edge-weight-dist-version.png as well as the ' + 
                  'default outputs/version/weighted/plots/edge-weights/edge-weight-dist-version.png.')
#parser.add_option('-o', '--out-file', type='string', default="viz/assays/zscore-rectfs-assays.png",
#        help='path/to/output_file.png. Default:')
parser.add_option('', '--pdf', action='store_true',
        help='Also store a pdf of the figures')
(opts, args) = parser.parse_args()


for i, version in enumerate(opts.version):
    print("")
    print("-"*30)
    t_settings.set_version(version)
    chemicals = sorted(utils.readItemList("%s/chemicals.txt" % (t_settings.INPUTSPREFIX)))

    interactome = t_settings.INTERACTOME

    print("Getting the weight, cost and direction of each edge from the interactome %s" % (interactome))
    lines = utils.readColumns(interactome, 1,2,3,4)
    edge_weights = {(u,v): float(w) for u,v,w,d in lines}
    edge_dir = {(u,v): True if d.lower() in ['true','t','dir','directed'] else False \
                for u,v,w,d in lines}

    # get the evidence from get_interaction_evidence.py
    #evidence_file = getev.getEvidenceFile(evidence_version, t_settings.DATADIR)
    #edge_dir = getev.getEdgeDir(edge_weights.keys(), evidence_file, split_family_nodes=False, add_ev_to_family_edges=False)
    # TODO try just the directed edges
    #weights = [edge_weights[e] for e in edge_weights if edge_dir[e] is True]
    weights = edge_weights.values()
    # for now, don't show the costs from the edges of weight 0
    costs = [-log(max(0.000001, w)) for w in weights if w != 0]
    #costs = [cost for cost in costs if cost < 4.5]
    # also incorporate the penalty into the plot
    penalty = log(t_settings.EDGE_PENALTIES.get(version, 1))
    costs = [cost+penalty for cost in costs]

    # if opts.ev_file is not None:
    #     evidence_file = opts.ev_file[i]
    # else:
    #     evidence_version = "2017_01pathlinker" if '2017_01' in t_settings.VERSION else "2016_05pathlinker"
    #     # get the evidence from get_interaction_evidence.py
    #     evidence_file = getev.getEvidenceFile(evidence_version, t_settings.DATADIR)
    # # now get the directions from the evidence file
    # print("Getting the edge direction of all edges in the interactome")
    # edge_dir = getev.getEdgeDir(edge_weights.keys(), evidence_file, split_family_nodes=False, add_ev_to_family_edges=False)

    print("Getting the weight each edge in the %d response networks" % (len(chemicals)))
    network_weights = []
    # first plot the distribution of edge weights in the response networks
    # get all of the response network edges and their edge weight in the interactome
    edgelinker_file = "%s/edgelinker/%%s-paths.txt" % (t_settings.RESULTSPREFIX) 
    for chemical in tqdm(chemicals):
        edges = t_utils.getEdges(paths_file=edgelinker_file % (chemical), max_k=200, ties=True)
        #tqdm.write("%d edges" % len(edges))
        for edge in edges:
            network_weights.append(edge_weights[edge])

    # now plot the distribution of edge weights
    #out_file_name = "response-network-weight-dist-dir-k500-%s.png" % (version)
    out_file_name = "response-network-weight-dist-k200-%s.png" % (version)
    out_dir = "%s/plots/edge-weights/" % (t_settings.RESULTSPREFIX)
    t_utils.checkDir(out_dir)
    out_file = "%s/%s" % (out_dir, out_file_name)
    if opts.compare_versions:
        out_dir_compare_versions = "viz/version_plots/edge-weights/"
        t_utils.checkDir(out_dir_compare_versions)
        out_file_compare_versions = "%s/%s" % (out_dir_compare_versions, out_file_name)

    print("Plotting response network edge weight histogram to %s" % (out_file))

    fig, ax = plt.subplots()
    print("%s of edges with a weight > %s" % (len([w for w in network_weights if w > 0.95]) / float(len(network_weights)), 0.95))

    ax.hist(network_weights, bins=30)

    plt.suptitle("Weights of %d response networks\n%s" % (len(chemicals), version))
    ax.set_xlabel("Response network edge weights")
    ax.set_ylabel("Frequency")

    plt.tight_layout()

    plt.savefig(out_file)
    if opts.pdf:
        plt.savefig(out_file.replace('.png', '.pdf'))
    plt.close()

    if opts.compare_versions:
        out_file2 = out_file_compare_versions
        print("Copying '%s' to '%s'" % (out_file, out_file2))
        shutil.copy(out_file, out_file2)
        if opts.pdf:
            print("Copying '%s' to '%s'" % (out_file.replace('.png','.pdf'), out_file2.replace('.png','.pdf')))
            shutil.copy(out_file.replace('.png','.pdf'), out_file2.replace('.png','.pdf'))

    # now plot the distribution of weights and costs in the interactome
    out_file_name = "interactome-weight-dist-%s.png" % (version)
    out_file = "%s/%s" % (out_dir, out_file_name)
    if opts.compare_versions:
        t_utils.checkDir(out_dir_compare_versions)
        out_file_compare_versions = "%s/%s" % (out_dir_compare_versions, out_file_name)

    print("Plotting interactome edge weight histogram to %s" % (out_file))
    fig, (ax1, ax2) = plt.subplots(nrows=2)

    # got the idea for the bins from here: https://stackoverflow.com/a/12176344/7483950
    #binwidth=1
    #ax1.hist(rec_zscores, bins=np.arange(min(rec_zscores), max(rec_zscores) + binwidth, binwidth))
    #ax2.hist(tf_zscores, bins=np.arange(-21, max(tf_zscores) + binwidth, binwidth))
    dir_weights = []
    undir_weights = []
    for edge, w in edge_weights.items():
        directed = edge_dir[edge]
        if directed:
            dir_weights.append(w)
        else:
            undir_weights.append(w)
    
    ax1.hist([undir_weights, dir_weights], bins=40, stacked=True, label=['undirected', 'directed'])
    ax1.legend()
    #ax1.hist(weights, bins=40)
    ax2.hist(costs, bins=40)

    # change the ticks for the edge cost histogram
    #start, end = ax2.get_xlim()
    #ax2.xaxis.set_ticks(np.arange(0, 4, 0.25))

    plt.suptitle("Histogram of weights and costs for interactome version\n%s" % (version))
    ax1.set_xlabel("Edge weights")
    ax2.set_xlabel("Edge costs")
    ax1.set_ylabel("Frequency")
    ax2.set_ylabel("Frequency")
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    plt.savefig(out_file)
    if opts.pdf:
        plt.savefig(out_file.replace('.png', '.pdf'))
    plt.close()

    if opts.compare_versions:
        out_file2 = out_file_compare_versions
        print("Copying '%s' to '%s'" % (out_file, out_file2))
        shutil.copy(out_file, out_file2)
        if opts.pdf:
            print("Copying '%s' to '%s'" % (out_file.replace('.png','.pdf'), out_file2.replace('.png','.pdf')))
            shutil.copy(out_file.replace('.png','.pdf'), out_file2.replace('.png','.pdf'))
