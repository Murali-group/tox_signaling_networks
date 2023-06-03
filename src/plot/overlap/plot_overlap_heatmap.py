# script to create a heatmap of the p-values of the overlap of two chemical subnetworks
""" Compute a p-value of the # of overlapping proteins between two chemicals by comparing the # of overlapping proteins between chemical A and 100 random networks of # rec and tfs of chemical B
"""

from optparse import OptionParser
import sys
import os
from tqdm import tqdm
import utilsPoirel as utils
sys.path.append("./src")
import toxcast_utils as t_utils

import matplotlib
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
#from matplotlib.ticker import FixedLocator, FixedFormatter
import pandas
import seaborn as sns
import numpy as np
import pdb
# turn off the annoying matplotlib font debug messages
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
#logging.getLogger('matplotlib.font_manager').disabled = True
#import plotly.plotly as py
#import plotly.graph_objs as go

CLUSTERING_METHODS=['average','ward','single','complete','weighted','centroid','median']


def parseArguments(args):
    usage = 'plot_overlap_heatmap.py [options]\n'
    parser = OptionParser(usage=usage)

    ## Common options
    parser.add_option('','--chemicals',type='string',
                      help='File containing chemicals to overlap. Required')
    parser.add_option('', '--mapping-file', type='string', metavar='STR',
                      help='File used to map to a different namespace. Network/edge IDs (uniprot ids) should be in the first column with the other namespace (gene name) in the second')
    parser.add_option('-o','--out-dir', action='store',
                      help='output to place the output files in')
    parser.add_option('','--paths-dir',type='string',metavar='STR',
                      help='Dir of Cyclinker results. Uses the ranked-edges file Required')
    parser.add_option('-k','--k-limit',action='store', type='int', default=200,
                      help='k limit for the overlap')
    parser.add_option('','--clustering-method', action='store', default='average',
                      help="Method for clustering the overlap matrix (using scipy.cluster.hierarchy.linkage). \n\tValid options are: %s" % (', '.join(CLUSTERING_METHODS)))
    parser.add_option('', '--pdf', action='store_true',
                      help='Also store a pdf of the figures')
    parser.add_option('','--forced',action='store_true',default=False,
                      help='Compute the overlap even if they will overwrite existing files.  Default is that algorithms are not run if output files exist.')

    # new options
    parser.add_option('-w','--write-overlap', action='store_true',
                      help='compute the overlap of the proteins in the chemical response networks')
    parser.add_option('','--split-family-nodes', action='store_true',default=False,
                      help='split family nodes when computing the overlap. Currently only works for writing the overlap')
    parser.add_option('','--plot-overlap', action='store_true',
                      help='compute the jaccard index of the proteins of each chemical response network pair and plot a clustered heatmap')
    parser.add_option('','--plot-hits-overlap', action='store_true',
                      help='compute the jaccard index of the hits of each chemical pair and plot a clustered heatmap')
    parser.add_option('','--interactome',type='string',metavar='STR',
                      help='interactome file used to limit the hits to those that are in the network. If not given, all hits will be used')
    parser.add_option('','--overlap-pvals', action='store',
                      help='Plot a heatmap of the pvals from the protein overlap clustering using the overlap pvals file created from overlap_sig.py')

    # parse the command line arguments
    (opts, args) = parser.parse_args()

    #if not opts.chemicals and not opts.sig_cutoff:
    if not opts.chemicals: 
        print("Must specify --chemicals for chemicals to use. Exiting")
        parser.print_help()
        sys.exit(1)
    if not opts.paths_dir:
        print("--paths-dir required. Exiting")
        parser.print_help()
        sys.exit(1)
    #if opts.plot_hits_overlap and not opts.interactome:
    #    print("ERROR: --interactome required to get the hits for --plot-hits-overlap")
    #    sys.exit(1)
    if opts.clustering_method not in CLUSTERING_METHODS:
        print("ERROR: --clustering-method '%s' invalid. Valid options are: %s" % (opts.clustering_method, ', '.join(CLUSTERING_METHODS)))
        sys.exit(1)

    return opts


def main(args):
    ## parse arguments
    opts = parseArguments(args)

    global uniprot_to_gene
    uniprot_to_gene = utils.readDict(opts.mapping_file, 1, 2)

    chemicals = utils.readItemList(opts.chemicals, col=1)

    print("Chemicals to overlap: ", len(chemicals))
    prefix = "" if opts.split_family_nodes is False else "-split"
    if opts.split_family_nodes is True:
        print("Splitting family nodes")

    # make sure the output directories are setup
    out_dir = "%s/%s%s/" % (opts.out_dir, opts.clustering_method, prefix)
    t_utils.checkDir(out_dir)
    print("Writing to out_dir: %s" % (out_dir))

    chem_prots, chem_edges, chem_rectfs, chem_paths_per_prot, chem_inter_paths_per_prot, chem_paths_per_edge = compute_overlap(opts, chemicals)

    if opts.write_overlap:
        # add a prefix of 'split-' if the split_family_nodes option is set
        prefix = "" if opts.split_family_nodes is False else "split-"
        write_num_networks_per_prot_overall(opts.out_dir + prefix, chem_paths_per_prot, chem_paths_per_edge)
        # also write the overlap of the intermediate proteins
        prefix += "intermediate-" 
        write_num_networks_per_prot_overall(opts.out_dir + prefix, chem_inter_paths_per_prot)

    # each of these options use the protein-overlap clustered heatmap
    if opts.plot_overlap or opts.plot_hits_overlap or opts.overlap_pvals:
        output_prefix = "%sprots-" % (out_dir)
        prot_df = compute_jaccard_and_cluster(chem_prots, output_prefix, method=opts.clustering_method, forced=opts.forced)
        prot_columns = prot_df.columns.tolist()

    if opts.plot_overlap:
        # clustering is already performed in compute_jaccard_and_cluster
        plot_jaccard_overlap(prot_df, output_prefix, cluster=False, pdf=opts.pdf, forced=opts.forced)

        output_prefix = "%srec-tfs-" % (out_dir)
        df = compute_jaccard_and_cluster(chem_rectfs, output_prefix, ordering=prot_columns, forced=opts.forced)
        plot_jaccard_overlap(df, output_prefix, cluster=False, pdf=opts.pdf, forced=opts.forced)

        output_prefix = "%sedges-" % (out_dir)
        df = compute_jaccard_and_cluster(chem_edges, output_prefix, ordering=prot_columns, forced=opts.forced)
        plot_jaccard_overlap(df, output_prefix, cluster=False, pdf=opts.pdf, forced=opts.forced)

        #output_prefix = "%spaths-per-prot-" % (out_dir)
        #df = compute_jaccard_and_cluster(chem_prot_num_paths, output_prefix, ordering=prot_columns, forced=opts.forced)
        #plot_jaccard_overlap(df, output_prefix, cluster=False, pdf=opts.pdf, forced=opts.forced)

    if opts.plot_hits_overlap:
        toxcast_data = t_utils.loadToxcastData(interactome_file=opts.interactome) 
        # read in the hit prots from the parsed files, needs the interactome 
        chem_hits = t_utils.getToxcastHits(
            toxcast_data=toxcast_data)
        # limit the chemicals to those that were passed in
        chem_hits = {c: h for c, h in chem_hits.items() if c in chemicals}
        hits_overlap(chem_hits, out_dir, ordering=prot_columns,
                     clustering_method=opts.clustering_method,
                     forced=opts.forced, pdf=opts.pdf)

    if opts.overlap_pvals:
        output_prefix = "%s" % (out_dir)
        plot_overlap_pval_heatmap(opts.overlap_pvals, output_prefix, ordering=prot_columns, forced=opts.forced)


def compute_overlap(opts, chemicals):
#def get_network_metrics(chemicals, k_limit)
    cyclinker_output = opts.paths_dir+'/%s-paths.txt'
    print("Reading paths for each chemical network from: %s" % (cyclinker_output))

    out_dir = "%s/chemicals/" % (opts.out_dir)
    # just write these once for the 'all' setting
    #if "all" in opts.out_dir and not os.path.isdir(out_dir):
    if "all" in opts.out_dir:
        print("Writing the num_paths_per_prot to %s" % (out_dir))
        t_utils.checkDir(out_dir)

    # dictionary containing each chemical's proteins. 
    chem_prots = {}
    chem_paths_per_prot = {}
    chem_inter_paths_per_prot = {} 
    chem_edges = {}
    chem_paths_per_edge = {}
    chem_rectfs = {}
    # count the number of times each protein is in each subnetwork, then rank them.
    for chemical in tqdm(chemicals):
        #proteins = t_utils.getProteins(paths=cyclinker_output % chemical, max_k=opts.k_limit)
        num_paths_per_prot, num_paths_per_edge, rec_tfs, num_inter_paths_per_prot = t_utils.get_path_metrics(cyclinker_output % (chemical), opts.k_limit,
                                                                                                             ties=True, split_family_nodes=opts.split_family_nodes)
        chem_prots[chemical] = num_paths_per_prot.keys()
        chem_paths_per_prot[chemical] = num_paths_per_prot
        chem_inter_paths_per_prot[chemical] = num_inter_paths_per_prot
        chem_edges[chemical] = num_paths_per_edge.keys()
        chem_paths_per_edge[chemical] = num_paths_per_edge
        chem_rectfs[chemical] = rec_tfs

        if 'all' in opts.out_dir:
            pref = "" if opts.split_family_nodes is False else "split-"
            # also write the counts for each chemical
            curr_out_dir = "%s/%s" % (out_dir, chemical)
            t_utils.checkDir(curr_out_dir)
            out_file = "%s/%s/%spaths-per-prot.txt" % (out_dir, chemical, pref)
            with open(out_file, 'w') as out:
                out.write("#prot\tgene\t# paths per prot\n")
                out.write('\n'.join(["%s\t%s\t%d" % (prot, uniprot_to_gene[prot], num_paths_per_prot[prot]) for prot in num_paths_per_prot]) + '\n')
            out_file = "%s/%s/%sintermediate-paths-per-prot.txt" % (out_dir, chemical, pref)
            with open(out_file, 'w') as out:
                out.write("#prot\tgene\t# paths per prot\n")
                out.write('\n'.join(["%s\t%s\t%d" % (prot, uniprot_to_gene[prot], num_inter_paths_per_prot[prot]) for prot in num_inter_paths_per_prot]) + '\n')

#    # Combine the set of proteins with the number of paths it appears in
#    # This will allow for the computation of the number of overlapping paths
#    chem_prot_num_paths = {}
#    for chemical in chemicals:
#        prots_nums = set()
#        for prot in chem_paths_per_prot[chemical]:
#            for i in range(chem_paths_per_prot[chemical][prot]):
#                prots_nums.add("%s-%d" % (prot, i))
#        chem_prot_num_paths[chemical] = prots_nums

    return chem_prots, chem_edges, chem_rectfs, chem_paths_per_prot, chem_inter_paths_per_prot, chem_paths_per_edge


def compute_jaccard_and_cluster(chem_prots, output_prefix, method="average", ordering=None, forced=False):
    #out_file = "%sjaccard-indeces.txt" % (output_prefix)
    out_file = "%sjaccard-clustered.txt" % (output_prefix)
    if forced or not os.path.isfile(out_file):
        print("Writing overlap to '%s'" % (out_file))
        # first compute the jaccard index of every chemical pair 
        # then write the overlap of the hits (rec and TF) and take some sort of ratio
        # dictionary of chemA-chemB pairs as the keys and the jaccard index as the overlap
        overlaps = {chemA: {} for chemA in chem_prots}
        for i, chemA in enumerate(chem_prots):
            for j, chemB in enumerate(chem_prots):
                # create a symmetrical matrix of the overlaps so that they can be clustered correctly
                if i < j:
                    chemA_prots = set(chem_prots[chemA])
                    chemB_prots = set(chem_prots[chemB])
                    jaccard_index = len(chemA_prots.intersection(chemB_prots)) / float(len(chemA_prots.union(chemB_prots))) 
                    overlaps[chemA][chemB] = jaccard_index
                    overlaps[chemB][chemA] = jaccard_index

        # convert to pandas dataframe for better visualization
        df = pandas.DataFrame(overlaps, dtype='float')
        df = df.fillna(0)

        # If the indeces were passed in, don't re-cluster the output. 
        # if the ordergin (ordering) is passed in to the function, use that instead
        if ordering is None:
            # perform the clustering here, and then write the clustered matrix to a CSV
            cg = sns.clustermap(df, method=method)
            # get the ordering of the chemicals after clustering
            column_indeces = cg.dendrogram_col.reordered_ind
            df = re_order_dataframe(df, column_indeces=column_indeces)
        else:
            df = re_order_dataframe(df, new_columns=ordering)
        # write the jaccard indeces
        df.round(2).to_csv(out_file)

    else:
        print("reading in overlap file '%s'. Use --forced to overwite it." % (out_file))
        df = pandas.read_csv(out_file, index_col=0)

    return df


def plot_jaccard_overlap(overlaps, output_prefix, method="average", cluster=False, pdf=False, forced=False):
    # plot the jaccard indeces as a clustered heatmap
    out_file = "%sjaccard-heatmap.png" % (output_prefix)

    if os.path.isfile(out_file) and not forced:
        print("%s already exists. Use --forced to overwrite it" % (out_file))
        return

    cmap = sns.cubehelix_palette(dark=0.2, reverse=True, as_cmap=True)
    if cluster:
        print("Plotting clustered overlap heatmap of the jaccard indeces to the file %s" % (out_file))
        #sns.heatmap(df).savefig(out_file)
        # taken from here: https://stanford.edu/~mwaskom/software/seaborn/generated/seaborn.clustermap.html
        cg = sns.clustermap(overlaps, method=method, cmap=cmap)
        # remove the tick labels
        cg.ax_heatmap.yaxis.set_ticklabels([])
        cg.ax_heatmap.xaxis.set_ticklabels([])
        # could rotate the labels, but they still overlap because they're too big
        #plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        #plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=1)
        cg.savefig(out_file)
        if pdf:
            print("Also generating pdf: %s " % out_file.replace('.png', '.pdf'))
            cg.savefig(out_file.replace('.png', '.pdf'))

    #if heatmap:
    else:
        # start a new figure to plot a non-clustered heatmap
        print("Plotting overlap heatmap of the jaccard indeces to the file %s" % (out_file))
        fig, ax = plt.subplots()
        sns.heatmap(overlaps, cmap=cmap)
        # TODO remove the tick labels
        #plt.yticks(rotation=0)
        #plt.xticks(rotation='vertical')
        plt.savefig(out_file)
        if pdf:
            print("Also generating pdf: %s " % (out_file.replace('.png', '.pdf')))
            plt.savefig(out_file.replace('.png', '.pdf'))


def plot_overlap_pval_heatmap(overlap_pvals_file, output_prefix, ordering=None, forced=False):
    """ Function to plot the significance of the overlap between two chemicals. 
    Keep the positions of the heatmap the same, but show the significance of each square's overlap
    """
    # First get the overlap pvals
    print("reading in overlap pvals from file '%s'." % (overlap_pvals_file))
    df = pandas.read_csv(overlap_pvals_file, index_col=0)

    # Now re-order the dataframe to match the overlap-clustered df
    df = re_order_dataframe(df, new_columns=ordering)
    out_file = "%soverlap-pvals-clustered.txt" % (output_prefix)
    if forced or not os.path.isfile(out_file):
        df.round(2).to_csv(out_file)

    # Now plot it
    #out_file = "%sjaccard-overlap-clustered-heatmap.png" % (output_prefix)
    #out_file = "%s-overlap-sig-heatmap.png" % (output_prefix)
    out_prefix = "%soverlap-pvals-" % (output_prefix)
    plot_jaccard_overlap(df, out_prefix, cluster=False, forced=forced)


def hits_overlap(
        chem_hits, out_dir, ordering=None,
        clustering_method='average', forced=False, pdf=False):

    # plot clustered on their own
    # The different clustering methods are explained on scipy: https://docs.scipy.org/doc/scipy-0.18.1/reference/cluster.hierarchy.html
    output_prefix = out_dir + "hits-" 
    df = compute_jaccard_and_cluster(
        chem_hits, output_prefix, method=clustering_method, forced=forced)
    plot_jaccard_overlap(df, output_prefix, cluster=False, pdf=pdf, forced=forced)

    # plot them with the response network protein clustering
    output_prefix = out_dir + "hits-prot-clust-" 
    #compute_jaccard_and_plot(chem_hits, output_prefix, ordering=ordering, forced=forced)
    df = re_order_dataframe(df, new_columns=ordering)
    plot_jaccard_overlap(df, output_prefix, cluster=False, pdf=pdf, forced=forced)


def re_order_dataframe(df, column_indeces=None, new_columns=None):
    #pdb.set_trace()
    if column_indeces is not None:
        cols = df.columns.tolist()
        new_columns = []
        for col_index in column_indeces:
            new_columns.append(cols[col_index])
        #pdb.set_trace()
    # for some reason, the index was the integer of the chemical ID. Map it to the string
    df.index = df.index.map(str)
    #df = df[new_cols]
    # the matrix is symmetrical so also reindex the rows
    df = df.reindex(new_columns, columns=new_columns)
    return df


def write_num_networks_per_prot_overall(output_prefix, chem_paths_per_prot, chem_paths_per_edge=None):
    # the number of networks each protein is in
    num_networks_per_prot_overall = {}
    # the number of paths each protein participates in
    num_paths_per_prot_overall = {}
    # the number of networks each edge participates in
    num_networks_per_edge_overall = {}
    # the number of paths each edge participates in
    num_paths_per_edge_overall = {}

    for chemical in chem_paths_per_prot:
        num_paths_per_prot = chem_paths_per_prot[chemical]
        # loop through the proteins and add to the overall count of each protein
        for protein in num_paths_per_prot:
            if protein not in num_networks_per_prot_overall:
                num_networks_per_prot_overall[protein] = 0 
            num_networks_per_prot_overall[protein] += 1 
            if protein not in num_paths_per_prot_overall:
                num_paths_per_prot_overall[protein] = []
            num_paths_per_prot_overall[protein].append(num_paths_per_prot[protein])

        if chem_paths_per_edge is not None:
            num_paths_per_edge = chem_paths_per_edge[chemical]
            # loop through the edges and add to the overall count of each edge
            for edge in num_paths_per_edge:
                if edge not in num_networks_per_edge_overall:
                    num_networks_per_edge_overall[edge] = 0 
                num_networks_per_edge_overall[edge] += 1 
                if edge not in num_paths_per_edge_overall:
                    num_paths_per_edge_overall[edge] = []
                num_paths_per_edge_overall[edge].append(num_paths_per_edge[edge])

    # now map the protein names to gene names
    #csbdb_interface = csbdb.Interface()
    ## map to entrez ID (GeneID in Uniprot)
    #uniprot_to_gene = csbdb_interface.map_ids(num_networks_per_prot_overall.keys(), 'UniprotKB', 'GeneName', '9606')
    #uniprot_to_gene = {prot:list(uniprot_to_gene[prot])[0] for prot in uniprot_to_gene}

    #if opts.chemical_type:
    #    out_file = "%s%soverlapping-proteins.txt" % (opts.output_prefix,'-'.join(opts.chemical_type))
    #else:
    out_file = "%soverlapping-proteins.txt" % (output_prefix)
    print("writing overlapping proteins to: %s" % (out_file))
    with open(out_file, 'w') as out:
        # write the header line
        out.write("#protein\tgene\t# of networks (out of %d)\ttotal # of paths protein is in\tavg # of paths protein is in\n" % len(chem_paths_per_prot))
        for protein in sorted(num_networks_per_prot_overall, key=num_networks_per_prot_overall.get, reverse=True):
            out.write("%s\t%s\t%d\t%d\t%0.2f\n" % (protein, uniprot_to_gene[protein], 
                num_networks_per_prot_overall[protein], sum(num_paths_per_prot_overall[protein]), np.average(num_paths_per_prot_overall[protein])))

    if chem_paths_per_edge is not None:
        out_file = "%soverlapping-edges.txt" % (output_prefix)
        print("writing overlapping edges to: %s" % (out_file))
        with open(out_file, 'w') as out:
            # write the header line
            out.write("#tail\thead\t# of networks (out of %d)\ttotal # of paths edge is in\tavg # of paths edge is in\n" % len(chem_paths_per_edge))
            for edge in sorted(num_networks_per_edge_overall, key=num_networks_per_edge_overall.get, reverse=True):
                out.write("%s\t%s\t%s\t%s\t%d\t%d\t%0.2f\n" % (edge[0], edge[1], uniprot_to_gene[edge[0]], uniprot_to_gene[edge[1]], 
                    num_networks_per_edge_overall[edge], sum(num_paths_per_edge_overall[edge]), np.average(num_paths_per_edge_overall[edge])))


if __name__ == "__main__":
    main(sys.argv)
