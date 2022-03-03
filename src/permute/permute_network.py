#! /usr/bin/env python

__author__ = "Jeff Law (jeffl@vt.edu)"

import networkx as nx
import random
from src.utils import file_utils as utils
from tqdm import tqdm
from optparse import OptionParser


def main():
    opts = parse_args()
    print("reading edges from %s" % (opts.interactome))
    G = nx.DiGraph(name="interactome")
    #physG = nx.Graph()
    if opts.test_run:
        G = build_test_graph() 
    else:
        for u,v,w in tqdm(utils.readColumns(opts.interactome, 1,2,3)):
            G.add_edge(u,v, weight=float(w))
            #if len(G) > 10000:
            #    break
    #print(nx.info(G))
    if opts.undirected:
        G = G.to_undirected()
        print('undirected ppi: %d nodes, %d edges' %(G.number_of_nodes(), G.number_of_edges()))
    else:
        print('directed ppi: %d nodes, %d edges' %(G.number_of_nodes(), G.number_of_edges()))
    print("Generating %d permuted (random) network(s), iterating over all of the edges %d time(s) for each random network" % (opts.num_random_networks, opts.num_iterations))
    for network_index in tqdm(range(opts.num_random_networks)):
        # swaps the edges in a copy of the original graph
        #permG = G.copy()
        permG = G
        permG = permute_network(permG, swap_phys_sig_sep=opts.swap_phys_sig_sep, split_weight=opts.split_by_weight, num_iterations=opts.num_iterations, max_tries=100)
        #print(nx.info(permG))

        out_file = '%s/permuted-network%d.txt' % (opts.out_dir, network_index+1)
        print("Writing %s" % (out_file))
        nx.write_weighted_edgelist(permG, out_file, comments='#', delimiter='\t')

    print("Done")


#@profile
def permute_network(permG, swap_phys_sig_sep=False, split_weight=None, num_iterations=1, max_tries=100, get_largest_cc=True, edge_lists=None):
    """
    """
    # speed up the process by mapping the node names to integers
    index = 1
    node2int = {}
    global int2node
    int2node = {}
    for node in permG.nodes():
        node2int[node] = index
        int2node[index] = node
        index += 1
    permG = nx.relabel_nodes(permG,node2int)

    print("swapping all of the edges of all of the nodes while preserving the in- and out-degree of each node ")
    num_swaps_performed = 0
    if swap_phys_sig_sep:
        if edge_lists is None:
            phys_edges, sig_edges = split_graph_phys_sig(permG)
        else:
            phys_edges, sig_edges = edge_lists
            phys_edges = [(node2int[u], node2int[v]) for u,v in phys_edges]
            sig_edges = [(node2int[u], node2int[v]) for u,v in sig_edges]
        print("swapping physical and signaling interactions separately with %d phys_edges and %d sig_edges" % (len(phys_edges), len(sig_edges)))

        # swap everything as undirected. The directed edges will be re-directed afterwards
        #permG = permG.to_undirected()

        directed = True
        if split_weight:
            print("and swapping edges above and below the given cutoff separately")
            phys_edges_high_weight, phys_edges_low_weight = split_graph_by_weight(permG, phys_edges, split_weight=split_weight)
            sig_edges_high_weight, sig_edges_low_weight = split_graph_by_weight(permG, sig_edges, split_weight=split_weight)
            edge_lists = [(sig_edges_high_weight, directed),
                    (sig_edges_low_weight, directed),
                    (phys_edges_high_weight, not directed),
                    (phys_edges_low_weight, not directed)]
        else:
            edge_lists = [(sig_edges, directed), (phys_edges, not directed)]

        for i in tqdm(range(num_iterations)):
            #print('-'*20)
            num_swaps_performed += swap_edges_random(permG, edge_lists=edge_lists, max_tries=max_tries)

        # direct the directed or signaling edges
        #permG = permG.to_directed()
        #
        # re-set the sig_edges list
        #if split_weight:
        #    sig_edges = sig_edges_high_weight + sig_edges_low_weight
        #
        #for u,v in sig_edges:
        #    # remove the opposite direction of the signaling edge
        #    permG.remove_edge(v, u)

        # This didn't give much of a speed up
#        for i in range(len(phys_edges)):
#            u,v = phys_edges[i]
#            # add the edge weights back
#            permG.edge[u][v]['weight'] = sig_weights[i] 
#
#        # direct the directed or signaling edges
#        permG = permG.to_directed()
#        for i in range(len(sig_edges)):
#            u,v = sig_edges[i]
#            # remove the opposite direction of the signaling edge
#            permG.remove_edge(v, u)
#            # add the edge weights back
#            permG.edge[u][v]['weight'] = sig_weights[i] 
    else:
        for i in tqdm(range(num_iterations)):
            num_swaps_performed += swap_edges_random(permG, max_tries=max_tries)
    print("performed a total of %d edge swaps" % (num_swaps_performed))

    if get_largest_cc:
        # Get largest connected component
        CC = sorted(nx.connected_components(permG.to_undirected()), key=len, reverse=True)
        print('-'*20)
        print("number of connected_components:", len(CC))
        largestCCnodes = set(CC[0])
        nodesToRemove = set(permG.nodes()).difference(largestCCnodes)
        permG.remove_nodes_from(nodesToRemove)
        print('\tdirected permuted ppi (largest weakly connected component): %d nodes, %d edges' %(permG.number_of_nodes(), permG.number_of_edges()))

    # put the original node names back
    permG = nx.relabel_nodes(permG,int2node)

    return permG


def split_graph_phys_sig(G):
    # pull out the signaling (directed) interactions from the physical (undirected) interactions
    print("separating the signaling (directed) interactions from the physical (undirected) interactions")
    phys_edges = []
    sig_edges = []
    for u,v in G.edges():
        # if the edge v -> u is also present, then add the edge to the undirected physG graph
        if G.has_edge(v,u):
            # if u is greater than v, then the v->u edge was already added
            # only need to add one direction of the undirected edges
            if u < v:
                phys_edges.append((u,v))
        # otherwise add the edge to the sigG graph
        else:
            sig_edges.append((u,v))
    return phys_edges, sig_edges


def split_graph_by_weight(G, edges=None, split_weight=0.6):
    """
    :param: edges - a list of edges to split. The network is required in order to get the edge weights.
    :param: split_weight - edges with weight above the split point will be put into the high_weight graph, edges with weight below will be in low_weight
    """
    # pull out the signaling (directed) interactions from the physical (undirected) interactions
    print("separating the low and high weight edges with a splitting point of %s" % (str(split_weight)))
#    high_weight = nx.DiGraph(name="%s_high_weight" % (G.name))
#    low_weight = nx.DiGraph(name="%s_low_weight" % (G.name))
#    for u,v in G.edges():
#        if G.edge[u][v]['weight'] > split_weight:
#            high_weight.add_edge(u,v, weight=G.edge[u][v]['weight'])
#        else:
#            low_weight.add_edge(u,v, weight=G.edge[u][v]['weight'])
    high_weight = []
    low_weight = []
    if edges is None:
        edges = G.edges()
    for u,v in edges:
        if G.adj[u][v]['weight'] > split_weight:
            high_weight.append((u,v))
        else:
            low_weight.append((u,v))
    return high_weight, low_weight


def swap_edges_random(G, edge_lists=[], max_tries=100, verbose=False):
    """
    :param: G: edges to swap 
    :param: masterG: the master graph which contains physical and signaling edges. 
    This master graph allows you to swap signaling edges amongst each other while not duplicating any of the signaling edges.
    The goal here is to retain the structure or relationship of the bi-directional physical edges
    """
    # If edge_lists are not specified, set the only subgraph as a pointer to the G graph 
    if len(edge_lists) == 0:
        # If there are no edge lists specified, then simply treat the network as directed as we don't have to worry about an undirected edge overwriting a directed edge 
        directed = True
        edge_lists = [(list(G.edges()), directed)]
    num_swaps_performed = 0
    for edge_list, directed in edge_lists:
        #print("permuting the edges of the graph %s" % (subgraph.name))
        #node_degrees_before = {n: (len(G.neighbors(n)), len(in_neighbors(G,n))) for n in nodes}
        num_edges = len(edge_list)
        #print("num nodes:", len(nodes), "num_tail_nodes:", len(tail_nodes))
        # randomly swap all of the edges of all of the nodes
        for i in range(num_edges):
            tries = 0
            while tries < max_tries:
                tries += 1
                rand_index1 = int(round(random.random() * (num_edges - 1)))
                rand_index2 = int(round(random.random() * (num_edges - 1)))
                #print(len(edge_list), rand_index1, rand_index2)
                #print(edge_list[0])
                #print(edge_list[rand_index2])
                u1, v1 = edge_list[rand_index1]
                u2, v2 = edge_list[rand_index2]
                if not directed:
                    # if this is an undirected edge, flip a coin to see which version of the undirected edge to swap
                    # meaning: u1-v1 is the same as v1-u1. Should we perform the swap u1-v2, v1-u2 - or - v1-v2, u1-u2?
                    if random.random() >= 0.5:
                        v1, u1 = u1, v1
                    if random.random() >= 0.5:
                        v2, u2 = u2, v2

                # don't allow self loops and make sure u->y or x->v don't already exist
                if u2 == u1 or u2 == v1 or v2 == u1 or v2 == v1 \
                   or G.has_edge(u2,v1) or G.has_edge(u1,v2):
                    continue
                # if this edge is undirected, make sure the reverse edge also isn't present or it will be overwritten
                if not directed and (G.has_edge(v1,u2) or G.has_edge(v2,u1)):
                    continue
                # break out if this loop. The edge swap passed the criteria
                break
            else:
                # if there was no eligible edge to swap with, skip this edge
                continue

            # perform the swap
            num_swaps_performed += 1
            #if verbose:
            #print("replacing %s->%s and %s->%s with %s->%s and %s->%s" % (u1,v1,u2,v2,u1,v2,u2,v1))

            # now replace u1->v1 and u2->v2 with u1->v2 and u2->v1
            # keep the weight of the old edge
            data1 = G.adj[u1][v1]
            data2 = G.adj[u2][v2]
            #G.add_edge(u1,v2, data1)
            #G.add_edge(u2,v1, data2)
            # networkx 2.1 update
            G.add_edge(u1,v2, **data1)
            G.add_edge(u2,v1, **data2)
            G.remove_edge(u1,v1)
            G.remove_edge(u2,v2)
            # also update the edge list
            edge_list[rand_index1] = (u1,v2)
            edge_list[rand_index2] = (u2,v1)
            # if this edge is undirected, also swap the reverse edge
            if not directed:
                G.add_edge(v2,u1, **data1)
                G.add_edge(v1,u2, **data2)
                G.remove_edge(v1, u1)
                G.remove_edge(v2, u2)

    #print("performed %d edge swaps" % (num_swaps_performed))
    return num_swaps_performed
#    node_degrees_after = {n: (len(G.neighbors(n)), len(in_neighbors(G,n))) for n in nodes}
#    #print(node_degrees_before)
#    #print(node_degrees_after)
#    for n in nodes:
#        in_degree, out_degree = node_degrees_before[n]
#        in_degree2, out_degree2 = node_degrees_after[n]
#        if in_degree != in_degree2 or out_degree != out_degree2:
#                print("ERROR: degree of %s changed from %d,%d (in,out) to %d,%d!" % (n, in_degree, out_degree, in_degree2, out_degree2))
#                sys.exit()


def build_test_graph():
    G = nx.DiGraph()
    # undirected edges
    G.add_edge(1,2, weight=0.5)
    G.add_edge(2,1, weight=0.5)
    G.add_edge(3,4, weight=0.5)
    G.add_edge(4,3, weight=0.5)
    G.add_edge(5,6, weight=0.7)
    G.add_edge(6,5, weight=0.7)
    G.add_edge(7,8, weight=0.7)
    G.add_edge(8,7, weight=0.7)
    # directed edges:
    G.add_edge(1,3, weight=0.5)
    G.add_edge(2,7, weight=0.5)
    G.add_edge(4,12, weight=0.5)
    G.add_edge(8,9, weight=0.7)
    G.add_edge(8,10, weight=0.7)
    G.add_edge(8,11, weight=0.7)
    G.add_edge(8,12, weight=0.7)
    G.add_edge(12,5, weight=0.7)

    return G


def parse_args():
    parser = OptionParser()

    parser.add_option('-i','--interactome',type='string',default="/home/jeffl/svnrepo/data/interactomes/human/2016_05/2016-06-02-human-ppi-mod-sig-5_100-weighted.txt",
            help="Path to interactome to use. \n\tDefault = '/home/jeffl/svnrepo/data/interactomes/human/2016_05/2016-06-02-human-ppi-mod-sig-5_100-weighted.txt'")
    parser.add_option('-o','--out_dir',type='string',default=".",
            help="Output directory to write the random networks to. \n\tDefault = '.'")   
    parser.add_option('','--num-random-networks',type='int',default='1',
            help="# of random networks to make. \n\tDefault = 1")   
    parser.add_option('','--num-iterations',type='int',default='1',
            help="# of times to iterate through all edges and swap them. \n\tDefault = 1")   
    parser.add_option('','--undirected',action='store_true',
            help="swap everything as undirected")   
    parser.add_option('','--swap-phys-sig-sep',action='store_true',
            help="swap the physical (undirected) and signaling (directed) edges seperately from each other")   
    parser.add_option('','--split-by-weight',type='float',
            help="swap the weights above and below the given weight separately")   
    parser.add_option('','--test-run',action='store_true',
            help="run on a test graph. Also print out the edges that are being replaced.")

    # parse the command line arguments
    (opts, args) = parser.parse_args()

    return opts

if __name__ == "__main__":
    main()
