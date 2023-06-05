#!/usr/bin/python

# code to post toxcast results to graphspace

from optparse import OptionParser
import sys
# turn off the many DEBUG messages
import logging
logging.basicConfig(level=logging.WARNING)
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

from src.utils import file_utils as utils
from src.graphspace import post_to_graphspace_base as gs_base
#import pdb
# used to get the evidence of the edges
#sys.path.append('/home/jeffl/src/python/Interactomes/Human')
from src.graphspace import get_interaction_evidence as getev
#try:
#import csbdb
#except ImportError:
#    sys.stderr.write("ImportError: This script imports csbdb to map to different namespaces. " +
#                     "You can download the CSBDB code using git:\n\tgit clone 'https://github.com/Murali-group/CSBDB.git'" +
#                     "\nThen add CSBDB to your PYTHONPATH variable.\n")
#    sys.exit(1)
try:
    # used for coloring the nodes according to hit/nonhit values
    sys.path.append("./src")
    import toxcast_utils as t_utils
except ImportError:
    pass

## GraphSpace variables
#GROUP_OWNER = 'jeffl@vt.edu'

#GROUP='ToxCast'
GROUP = ''

# Filtered detection method (from parse-human-ppi.py)
METHODS_TO_IGNORE = ('MI:0036',  # domain fusion
                     'MI:0046',  # experimental knowledge based
                     'MI:0063',  # interaction prediction
                     'MI:0064',  # interologs mapping
                     'MI:0085',  # phylogenetic profile
                     'MI:0087',  # predictive text mining
                     'MI:0105',  # structure based prediction
                     'MI:0363',  # inferred by author
                     'MI:0364',  # inferred by curator
                     'MI:0686',  # unspecified method coexpression
                     'MI:0045',  # experimental interaction detection
                     'MI:0813',  #proximity litigaion assay:Aditya
)

# global variables will be populated by main()
DATADIR = ''
PPI = ''
PPIEDGES = []
PPIWEIGHTS = {}

# Dictionaries of node and edge properties
NODE_COLORS = {
        'target' : '#FFFF60',  # yellow
        'source' : '#8CE1DC',  # blue
        'default' : '#D8D8D8',  # gray
        'intermediate_rec' : '#D8D8D8',  # gray
        'intermediate_tf' : '#D8D8D8',  # gray
        # for the parent/compound nodes
        'Targets' : '#FFFF60',  # yellow
        'Sources' : '#8CE1DC',  # blue
        'Intermediate Proteins' : '#D8D8D8',  # gray
        'Intermediate Receptors' : '#D8D8D8',  # gray
        'Intermediate TFs' : '#D8D8D8',  # gray
         #'netpath_kegg_egfr' : '#AC58FA',  # purple
         #'responsive' :  "#FACC2E",  # orange
         #'non-responsive' :  "#AC58FA",  # purple
}
CASESTUDY_NODE_COLORS = {
        'target' : '#D8D8D8',  # gray
        'source' : '#D8D8D8',  # gray
}

NODE_SHAPES = {
        'source'           : 'diamond',
        'target'           : 'rectangle',
        'default'          : 'ellipse',
        'intermediate_rec' : 'triangle',
        'intermediate_tf'  : 'rectangle',
}

EDGE_DIR_MAP = {
        'physical':False,
        'phosphorylation':True,
        'enzymatic':True,
        'activation':True,
        'inhibition':True,
        'spike_regulation':True,
        'dir_predicted': True,
}

EDGE_COLORS = {
        'physical': '#27AF47',  # green
        'phosphorylation': '#F07406',  # orange
        'enzymatic': 'brown',  # blue
        #'enzymatic': '#DD4B4B',  # red
        'activation': 'grey',
        'inhibition': 'grey',
        'spike_regulation': 'red',
        'dir_predicted': '#2A69DC',
}

CASESTUDY_EDGE_COLORS = {
        'physical': '#27AF47',  # green
        'phosphorylation': '#F07406',  # orange
        'enzymatic': 'brown',  # blue
        'activation': 'grey',  # orange
        'inhibition': 'grey',  # orange
        'spike_regulation': 'grey',  # orange
        'dir_predicted': '#2A69DC',
}
# order to use when deciding which color an edge will get
edge_type_order = ['phosphorylation', 'enzymatic', 'spike_regulation', 'activation', 'inhibition', 'dir_predicted', 'physical']


# CSBDB interface for fetching evidence
#db = csbdb.Interface()
# use the most recent proteome (2017_01)
#db.current_schema_name = 'build_reviewed_proteome_2016_05'
#db.current_schema_name = 'build_full_proteome_2016_05'

# mapping from uniprot ID gene name (including family nodes)
uniprot_to_gene = {}


def main(args):
    global PPI, PPIEDGES, PPIWEIGHTS
    global EVIDENCEFILES, DATADIR
    global CASESTUDY, PARENTNODES

    opts, args = parseArgs(args)

    # (0) Determine map file, PPI, and versions
    DATADIR = opts.datadir
    CASESTUDY = opts.casestudy
    if CASESTUDY is True:
        print("Using 'case study' colors and borders")
    PARENTNODES = opts.include_parent_nodes

    PPI = opts.ppi
    lines = utils.readColumns(PPI,1,2,3)
    PPIEDGES = [(u,v) for u,v,w in lines]
    PPIWEIGHTS = {(u,v):float(w) for u,v,w in lines}

    if opts.paths or opts.ranked_edges:
        prededges = readNetwork(paths=opts.paths, ranked_edges=opts.ranked_edges, k_limit=opts.k_limit)

    lines = utils.readColumns(opts.sourcetarget,1,2)
    sources = set([acc for acc, node_type in lines if node_type.lower() in ['source', 'receptor']])
    targets = set([acc for acc, node_type in lines if node_type.lower() in ['target', 'tf']])
    human_rec = set() 
    human_tfs = set()
    if opts.human_rec_tfs is not None:
        # also get the human rec and tfs
        lines = utils.readColumns(opts.human_rec_tfs,1,2)
        human_rec = set([acc for acc, node_type in lines if node_type.lower() in ['source', 'receptor']])
        human_tfs = set([acc for acc, node_type in lines if node_type.lower() in ['target', 'tf']])

    global sig_chemicals
    sig_chemicals = utils.readItemSet(opts.sig_chemicals)
    print("Posting the nodes/edges that are common to %d of the %d significant networks" % (int(opts.overlap_cutoff*len(sig_chemicals)), len(sig_chemicals)))

    #if opts.chemicalID:
    #global CHEMICALS
    #CHEMICALS = t_utils.loadChemicalMap(PPI)
    #chemDSStoName, chemNametoDSS = t_utils.getChemicalNameMaps()

    # get attributes of nodes and edges from the graph_attr file
    graph_attr = {}
    # description of a node/edge style attribute that will be added to the node/edge popup 
    # two-level dictionary structure: 
    # node/edge
    #   attribute
    #     description
    attr_desc = {}
    if opts.graph_attr:
        # keep the order of the pathways by order of highest posterior probability
        #pathway_colors = collections.OrderedDict()
        print("Adding graph attributes from '%s' " % (opts.graph_attr))
        # example node attribute:
        # color blue    p1|p2|p3
        # example edge attribute:
        # edge_style dotted    p1-p2|p2-p3
        # example compound node. Here p1, p2, and p3 will have the parent attribute set to 'parent1' (i.e. they will belong to the same compound node parent1)
        # parent    parent1  p1|p2|p3
        # then to set the attributes of 'parent1', specify it as the node
        # color blue    parent1
        for style, style_attr, prots, desc in utils.readColumns(opts.graph_attr, 1,2,3,4):
            for prot in prots.split('|'):
                if prot not in graph_attr:
                    graph_attr[prot] = {}
                graph_attr[prot][style] = style_attr
            attr_desc[(style, style_attr)] = desc
            #graph_attributes[group_number] = {"style": style, "style_attr": style_attr, "prots": prots.split(','), "desc":desc}

    # add the edge overlap as a graph attribute
    if opts.edge_overlap is not None:
        num_paths_per_edge = {}
        num_nets_per_edge = {}
        for u, v, num_networks, num_paths in utils.readColumns(opts.edge_overlap, 1, 2, 5, 6):
            if float(num_networks) / float(len(sig_chemicals)) <= opts.overlap_cutoff:
                break
            num_paths_per_edge[(u,v)] = int(num_paths) 
            num_nets_per_edge[(u,v)] = int(num_networks) 

        print("Total number of edges: %s" % (len(num_paths_per_edge)))
        # use these edges as the "predicted" edges
        prededges = {}
        # use the number of paths as the k filter
        k = 1
        # sort by the number of paths in decreasing order
        for edge in sorted(num_paths_per_edge, key=num_paths_per_edge.get, reverse=True):
            prededges[edge] = k
            k += 1

        #normalized_num_paths_per_edge = {}
        # now normalize the edge width between 2 and 15 (seemed to be the best sizes)
        a = 2
        #b = 15
        b = 12
        max_num_paths = max(num_paths_per_edge.values())
        min_num_paths = min(num_paths_per_edge.values())
        print("edge max_num_paths = %s, min_num_paths = %s" % (max_num_paths, min_num_paths))
        for u,v in num_paths_per_edge:
            #normalized_num_paths_per_edge[(u,v)] = (b-a) * ((num_paths_per_edge[(u,v)] - min_num_paths) / (max_num_paths - min_num_paths)) + a
            normalized_num_paths = (b-a) * (float(num_paths_per_edge[(u,v)] - min_num_paths) / float(max_num_paths - min_num_paths)) + a
            edge_str = "%s-%s" % (u,v)
            if edge_str not in graph_attr:
                graph_attr[edge_str] = {}
                attr_desc[edge_str] = {}
            graph_attr[edge_str]['width'] = normalized_num_paths
            attr_desc[edge_str]["# paths edge is in"] = num_paths_per_edge[(u,v)]
            attr_desc[edge_str]["# networks edge is in"] = num_nets_per_edge[(u,v)]

    # also add the node overlap as a node attribute
    if opts.node_overlap is not None:
        num_paths_per_node = {}
        num_nets_per_node = {}
        for n, num_networks, num_paths in utils.readColumns(opts.node_overlap, 1, 3, 4):
            #if int(num_networks) <= opts.min_num_networks:
            if float(num_networks) / float(len(sig_chemicals)) <= opts.overlap_cutoff:
                break
            num_paths_per_node[n] = int(num_paths) 
            num_nets_per_node[n] = int(num_networks) 

        print("Total number of nodes: %s" % (len(num_paths_per_node)))
        # now normalize the width and height between 20 and 150 (seemed to be the best)
        a = 20
        #b = 150
        b = 120
        max_num_paths = max(num_paths_per_node.values())
        min_num_paths = min(num_paths_per_node.values())
        print("node max_num_paths = %s, min_num_paths = %s" % (max_num_paths, min_num_paths) )
        for n in num_paths_per_node:
            #normalized_num_paths_per_node[e] = (b-a) * ((num_paths_per_node[e] - min_num_paths) / (max_num_paths - min_num_paths)) + a
            normalized_num_paths = (b-a) * (float(num_paths_per_node[n] - min_num_paths) / float(max_num_paths - min_num_paths)) + a
        # max_num_nets = max(num_nets_per_node.values())
        # min_num_nets = min(num_nets_per_node.values())
        # print("node max_num_nets = %s, min_num_nets = %s" % (max_num_nets, min_num_nets) )
        # for n in num_nets_per_node:
        #     #normalized_num_nets_per_node[e] = (b-a) * ((num_nets_per_node[e] - min_num_nets) / (max_num_nets - min_num_nets)) + a
        #     normalized_num_nets = (b-a) * (float(num_nets_per_node[n] - min_num_nets) / float(max_num_nets - min_num_nets)) + a
            if n not in graph_attr:
                graph_attr[n] = {}
                attr_desc[n] = {}
            graph_attr[n]['width'] = normalized_num_paths
            graph_attr[n]['height'] = normalized_num_paths
            attr_desc[n]["# paths node is in"] = num_paths_per_node[n]
            attr_desc[n]["# networks node is in"] = num_nets_per_node[n]
            num_hits = 0
            # TODO if this is a hit, also write the # of networks in which the node is hit 
            #for chem in sig_chemicals:
            #    for split_n in n.split(','):
            #        if split_n in CHEMICALS.chemical_protein_hit[chem]:
            #            num_hits += 1
                #if n == "P00533" and n not in CHEMICALS.chemical_protein_hit[chem]:
                #    print("%s (%s) has EGFR as an intermediate node" % (chemDSStoName[chem], chem))
            if num_hits > 0:
                attr_desc[n]["# hits"] = num_hits

    #sys.exit()
    global uniprot_to_gene
    #uniprot_to_gene = utils.readDict(getev.getMappingFile(opts.version, opts.datadir), 1, 2)
    #gene_to_uniprot = utils.readDict(getev.getMappingFile(opts.version, opts.datadir), 2, 1)
    #uniprot_to_gene = {} 
    if opts.mapping_file is not None:
        uniprot_to_gene = utils.readDict(opts.mapping_file, 1, 2)

    if opts.ctd_support is not None:
        print("Getting CTD support counts from %s" % (opts.ctd_support))
        num_intxs_per_gene = {}
        num_intxs_per_node = {}
        for cas, gene, interaction_action, pubmedids in utils.readColumns(opts.ctd_support, 3, 4, 10, 11):
            if 'phosphorylation' in interaction_action:
                if gene not in num_intxs_per_gene:
                    num_intxs_per_gene[gene] = 0
                num_intxs_per_gene[gene] += 1 

        prednodes = set([N for edge in prededges for N in edge])

        for N in prednodes:
            for n in N.split(","):
                gene = uniprot_to_gene[n]
                if gene in num_intxs_per_gene:
                    # for now, take the node in the family node with the maximum support
                    if N in num_intxs_per_node and num_intxs_per_node[N] > num_intxs_per_gene[gene]:
                        pass
                    else:
                        num_intxs_per_node[N] = num_intxs_per_gene[gene]
                    #print(uniprot_to_gene[N], gene, num_intxs_per_gene[gene] )
        print("%d out of %d nodes have support in CTD" % (len(num_intxs_per_node), len(prednodes)))

        # write these counts to a table
        out_file = "%s-ctd-support.txt" % (opts.node_overlap)
        print("Writing support counts to %s" % (out_file))
        with open(out_file, 'w') as out:
            out.write("#uniprot\tgene\tmax_num_phospho_intxs\n")
            # write the sorted results to a file
            out.write('\n'.join(["%s\t%s\t%d" % (N, uniprot_to_gene[N], num_intxs_per_node[N]) for N in sorted(num_intxs_per_node, key=num_intxs_per_node.get, reverse=True)]) + '\n')

        # most of the nodes have < 500, but the range is up to 5000 
        max_support = 100
        #min_support = 30
        b = 1
        a = 0
        for N in prednodes:
            if N in num_intxs_per_node:
                num_intxs = num_intxs_per_node[N]
                if num_intxs > max_support:
                    num_intxs = max_support
                normalized_support = (b-a) * (float(num_intxs - 0) / float(max_support - 0)) + a
            else:
                normalized_support = 0
                num_intxs = 0 
            graph_attr[N]['background-opacity'] = normalized_support
            # try using this color as the darkest
            graph_attr[N]['background-color'] = '#4286f4'
            graph_attr[N]['text-valign'] = 'bottom'
            attr_desc[N]["# CTD phosphorylation interactions"] = num_intxs_per_node[N] if N in num_intxs_per_node else 0

    # set the case study colors
    if CASESTUDY is True:
        if opts.graph_attr is not None:
            # if there are other colors present, then make the sources and targets gray by default because they could have other colors
            NODE_COLORS.update(CASESTUDY_NODE_COLORS)
        EDGE_COLORS.update(CASESTUDY_EDGE_COLORS)

    # Now post to graphspace!
    G = constructGraph(opts, prededges, sources, targets, extra_rec=human_rec, extra_tfs=human_tfs, graph_attr=graph_attr, attr_desc=attr_desc)

    # TODO fix mapping of family nodes
    #G = relabelNodes(G)
    # rename the nodes to have Gene names rather than Uniprot IDs
    #G = nx.relabel_nodes(G, mapNodes(G.nodes(), from_namespace='uniprotkb', to_namespace='genename'))

    # post to graphspace
    gs = GraphSpace(opts.username, opts.password)
    gs_graph = gs.get_graph(opts.graph_name, owner_email=opts.username)
#    if gs_graph is None:
    if opts.apply_layout is None:
        print("\nPosting graph '%s' to graphspace\n" % (opts.graph_name))
        gs_graph = gs.post_graph(G)
    else:
        # Currently layouts store style attributes as well as x and y positions
        # so if you layed out the graph how you wanted it but now want to update the style (such as edge width), 
        # you will need to copy the x and y positions of the layout you made to the updated graph
        # I created the layout 'layout1', so I set that as the default
        layout_name = opts.layout_name
        if opts.apply_layout is not None:
            # if a layout was created for a different graph name, try to copy that layout here
            print("checking if layout '%s' exists for graph %s" % (layout_name, opts.apply_layout))
            layout = gs.get_graph_layout(graph_name=opts.apply_layout,layout_name=layout_name)
        else:
            print("checking if layout '%s' exists for this graph (%s)" % (layout_name, opts.graph_name))
            layout = gs.get_graph_layout(graph=gs_graph,layout_name=layout_name)
        if layout is not None:
            # set the x and y position of each node in the updated graph to the x and y positions of the layout you created
            print("Setting the x and y coordinates of each node to the positions in %s" % layout_name)
            for node, positions in layout.positions_json.items():
                G.set_node_position(node_name=node, x=positions['x'], y=positions['y'])

        # TODO fix this
        print("\nPosting graph '%s' to graphspace\n" % (opts.graph_name))
        gs_graph = gs.post_graph(G)
        # # now re-post the graph and the positions should be set
        # print("\nGraph '%s' already exists. Updating it\n" % (opts.graph_name))
        # gs_graph = gs.update_graph(G, graph_name=opts.graph_name, owner_email=USERNAME)
    print(gs_graph.url)

    if opts.group is not None:
        # share with group 'ToxCast'    
        # you can also share your graph with a group
        # create the group if it doesn't exist
        #group = gs.post_group(GSGroup(name='icsb2017', description='sample group'))
        # or get the group you already created it
        print("sharing graph with group '%s'" % (opts.group))
        group = gs.get_group(group_name=opts.group)
        #print(group.url)
        gs.share_graph(graph=gs_graph, group=group)


def readNetwork(paths=None, ranked_edges=None, k_limit=200):
    """ Read the PathLinker paths or ranked_edges output. 
        Get all of the edges that have a k less than the k_limit.
    """
    if paths:
        # Predicted paths from pathlinker
        lines = utils.readColumns(paths,1,2,3)
        prededges = {}
        edges = set()
        for k, path_score, path in lines:
            # get all of the edges in the paths that have a k less than the k_limit
            if int(k) > k_limit:
                break

            path = path.split('|')

            for i in range(len(path)-1):
                edge = (path[i], path[i+1])
                if edge not in edges:
                    edges.add(edge)
                    prededges[edge] = int(k)

    if ranked_edges:
        # Predicted edges from pathlinker
        lines = utils.readColumns(ranked_edges,1,2,3)
        # get all of the edges that have a k less than the k_limit
        prededges = {(u,v):int(k) for u,v,k in lines if int(k) <= k_limit}

    return prededges


def constructGraph(opts, prededges, sources, targets, extra_rec=[], extra_tfs=[], graph_attr={}, attr_desc={}):
    '''
    Posts the toxcast pathlinker result to graphspace

    :param source: list of source nodes
    :param targets: list of target nodes
    :param graphid: name of graph to post
    :param outfile: output JSON file that will be written
    '''
    # NetworkX object
    #G = nx.DiGraph(directed=True)
    G = GSGraph()
    edges_file = opts.ranked_edges if opts.ranked_edges is not None else opts.paths
    # get metadata
    desc = getGraphDescription(edges_file, opts.ppi, pathway_colors=None)
    metadata = {'description':desc,'tags':[], 'title':''}
    if opts.tag:
        metadata['tags'] = opts.tag
    #if opts.chemicalID:
    #    metadata['title'] = CHEMICALS.chemDSStoName[opts.chemicalID]

    evidence_file = opts.ev_file
    if evidence_file is None:
        # get the evidence from get_interaction_evidence.py
        evidence_file = getev.getEvidenceFile('2018_01pathlinker', opts.datadir)
    # TODO make it an option to add evidence to family edges
    add_ev_to_family_edges = True 
    family_ppi_evidence = None 
    if add_ev_to_family_edges:
        evidence, family_ppi_evidence, edge_types, edge_dir = getev.getEvidence(prededges.keys(), evidence_file=evidence_file, split_family_nodes=True, add_ev_to_family_edges=add_ev_to_family_edges)
    else:
        evidence, edge_types, edge_dir = getev.getEvidence(prededges.keys(), evidence_file=evidence_file, add_ev_to_family_edges=add_ev_to_family_edges)

    prednodes = set([t for t,h in prededges]).union(set([h for t,h in prededges]))

    # set of parent nodes to add to the graph
    parents_to_add = {}

    ## add GraphSpace/Cytoscape.js attributes to all nodes.
    for n in prednodes:
        #default is gray circle
        node_type = 'default'
        # default parent
        parent = 'Intermediate Proteins'
        if n in sources:
            # if n is the source, make it a blue triangle
            node_type = 'source'
            # set the parent node for this source 
            parent = 'Sources'
        elif n in targets:
            # if n is a taret, make it a yellow square
            node_type = 'target'
            parent = 'Targets'
        elif n in extra_rec:
            # gray triangle
            node_type = 'intermediate_rec'
            parent = 'Intermediate Receptors'
        elif n in extra_tfs:
            # gray square
            node_type = 'intermediate_tf'
            parent = 'Intermediate TFs'

        if n not in graph_attr:
            graph_attr[n] = {}
            graph_attr[n]['parent'] = parent
        # only overwrite the parent if it is a source or target
        elif node_type == 'source' or node_type == 'target':
            graph_attr[n]['parent'] = parent

        edgeswithnode = set([(t,h) for t,h in prededges if t==n or h==n])
        pathswithnode = set([int(prededges[e]) for e in edgeswithnode])
        k_value = min(pathswithnode)
        # set the name of the node to be the gene name and add the k to the label
        gene_name = uniprot_to_gene[n]
        # make sure the gene name is sorted
        gene_name = ','.join(sorted(gene_name.split(',')))
        short_name = gene_name
        # TODO if the the family name is too long, then just choose one of the genes and add -family to it
        if len(gene_name) > 10:
            short_name = "%s-family" % (gene_name.split(',')[0])

        attr_dict = {}
        if PARENTNODES is True:
            # set the parent if specified
            if n in graph_attr and 'parent' in graph_attr[n]:
                # set the parent of this node
                parent = graph_attr[n]['parent']
                attr_dict['parent'] = parent
            # also add this parent to the set of parent/compound nodes to add to the graph
            # keep track of the lowest k value so it will work correctly with the sliding bar
            if parent not in parents_to_add or k_value < parents_to_add[parent]:
                parents_to_add[parent] = k_value

        #node_popup = getNodeAnnotation(n, pathswithnode, attr_desc=attr_desc)
        node_popup = gs_base.buildNodePopup(n, attr_val=attr_desc)

        label = "%s"%(short_name)

        # TODO set the node label smaller than the gene name for large family node labels
        G.add_node(gene_name, attr_dict=attr_dict, popup=node_popup,
                   label=label, k=k_value)

        attr_dict = {}
        shape = NODE_SHAPES[node_type]
        color = NODE_COLORS[node_type]
        style = 'solid'
        width = 45
        height = 45
        border_width = 2
        border_color = "#7f8184"
        bubble = NODE_COLORS[node_type]
        graph_attr[n]['text-valign'] = 'bottom'
        if n in graph_attr:
            # allow for any attribute to be set
            for s in graph_attr[n]:
                attr_dict[s] = graph_attr[n][s]
            if 'color' in graph_attr[n]:
                color = graph_attr[n]['color']
            if 'shape' in graph_attr[n]:
                shape = graph_attr[n]['shape']
            if 'style' in graph_attr[n]:
                style = graph_attr[n]['style']
            if 'border_color' in graph_attr[n]:
                border_color= graph_attr[n]['border_color']
            if 'border_width' in graph_attr[n]:
                border_width= graph_attr[n]['border_width']
            if 'width' in graph_attr[n]:
                width = graph_attr[n]['width']
            if 'height' in graph_attr[n]:
                height = graph_attr[n]['height']
            if 'bubble' in graph_attr[n]:
                bubble = graph_attr[n]['bubble']
        #border_color = color if border_color is None else border_color
        #bubble = color if bubble is None else bubble

        #if CASESTUDY is False:
            # # add the double borders if this isn't a casestudy
            # #if opts.chemicalID and n in CHEMICALS.chemical_protein_hit[opts.chemicalID] \
            # #        and (CASESTUDY is False or (node_type == "source" or node_type == "target")):
            # # if this is tested for any chemical, but a double border around it
            # # TODO add the actual counts to the nodes
            # for chem in sig_chemicals:
            #     if border_color == "maroon":
            #         break
            #     for split_n in n.split(','):
            #         if split_n in CHEMICALS.chemical_protein_hit[chem]:
            #     #print(n, CHEMICALS.chemicals[chemicalID]['accs'][n])
            #             style = 'double'
            #             border_width = 10
            #             if CHEMICALS.chemical_protein_hit[chem][split_n] > 0:
            #                 border_color = 'maroon'
            #             else: 
            #                 border_color = 'black'
            # if this is a source or target, then apply the double border automatically
            #if opts.chemicalID and (node_type in ["source", "target"]): 
            # if node_type in ["source", "target"]: 
            #     style = 'double'
            #     border_width = 10
            #     border_color = 'maroon'

        G.add_node_style(gene_name, shape=shape, attr_dict=attr_dict, color=color, width=width, height=height,
                         style=style, border_color=border_color, border_width=border_width, bubble=bubble)

    # now add the parent nodes
    for parent, k_value in parents_to_add.items():
        # TODO add a popup for the parent which would be a description of the function or pathway
        if ('parent', parent) in attr_desc:
            popup = attr_desc[('parent', parent)]
        else:
            popup = parent

        if parent in graph_attr and 'label' in graph_attr[parent]:
            label = graph_attr[parent]['label']
        else:
            label = parent.replace("_", " ")

        #popup = parent
        # leave off the label for now because it is written under the edges and is difficult to see in many cases
        # I requested this feature on the cytoscapejs github repo: https://github.com/cytoscape/cytoscape.js/issues/1900
        G.add_node(parent, popup=popup, k=k_value, label=label)

        # the only style we need to set for the parent/compound nodes is the color
        # default color for sources and targets, and intermediate nodes
        if parent in NODE_COLORS:
            color = NODE_COLORS[parent]
        if parent in graph_attr and 'color' in graph_attr[parent]:
            color = graph_attr[parent]['color']
        attr_dict = {} 
        # set the background opacity so the background is not the same color as the nodes
        attr_dict["background-opacity"] = 0.3
        attr_dict['font-weight'] = "bolder"
        attr_dict['font-size'] = "24px"
        attr_dict['text-outline-width'] = 3
        # remove the border around the parent nodes 
        border_width = 0
        valign="bottom"

        G.add_node_style(parent, attr_dict=attr_dict, color=color, border_width=border_width, bubble=color, valign=valign)

    # Add all of the edges and their Graphspace/Cytoscape.js attributes
    for (u,v) in prededges:
        # get the main edge type
        main_edge_type = getMainEdgeType(edge_types[(u,v)])

        # add edge k
        k_value = prededges[(u,v)]
        if 'activation' not in edge_types[(u,v)] and 'inhibition' in edge_types[(u,v)]:
            # I don't think the graphspace interface has anything for this, so add it here
            arrow_shape = 'tee'
        else:
            arrow_shape = 'triangle'

        gene_name_u = ','.join(sorted(uniprot_to_gene[u].split(',')))
        gene_name_v = ','.join(sorted(uniprot_to_gene[v].split(',')))

        # family_ppi_evidence will be None if we are not including famliy edges
        edge_popup = getEdgeAnnotation(u,v,k_value, evidence, family_ppi_evidence=family_ppi_evidence, attr_desc=attr_desc)
        G.add_edge(gene_name_u,gene_name_v,directed=edge_dir[(u,v)],popup=edge_popup,k=k_value)

        attr_dict = {}
        attr_dict['opacity'] = 0.7
        color = EDGE_COLORS[main_edge_type]
        edge_style = 'solid'
        width = 1.5
        edge_str = "%s-%s" % (u,v)
        if edge_str in graph_attr:
            if 'color' in graph_attr[edge_str]:
                color = graph_attr[edge_str]['color']
            if 'arrow_shape' in graph_attr[edge_str]:
                arrow_shape= graph_attr[edge_str]['arrow_shape']
            if 'edge_style' in graph_attr[edge_str]:
                edge_style = graph_attr[edge_str]['edge_style']
            if 'width' in graph_attr[edge_str]:
                width = graph_attr[edge_str]['width']

        G.add_edge_style(gene_name_u, gene_name_v, attr_dict=attr_dict,
                         directed=EDGE_DIR_MAP[main_edge_type], color=color, width=width, arrow_shape=arrow_shape,edge_style=edge_style)
        G.set_data(metadata)
        G.set_name(opts.graph_name)
    return G


def getMainEdgeType(edge_types):
    """ a single edge can have multiple edge types according to the different sources or databases
    Choose a main edge type here
    *edge_types* the set of edge types for a given edge 
    """
    main_edge_type = None
    for edge_type in edge_type_order:
        if edge_type in edge_types:
            main_edge_type = edge_type 
            break

    if main_edge_type is None:
        sys.stderr.write("WARNING: %s has no edge type. edge_types[%s]: %s. " % (str((u,v)),str((u,v)), str(edge_types[(u,v)])))
        sys.stderr.write("Evidence for edge: %s\n" % (str(evidence[(u,v)])))
        sys.stderr.write("\tSetting to 'activation'\n")
        #CHECK
        #sys.exit(1)
        # set this as the default for now
        main_edge_type = "activation"

    #if main_edge_type == '':
    #    raise NoEdgeTypeError("ERROR: %s,%s has no edge type. edge_types[%s,%s]: %s" % (u,v,u,v, str(edge_types[(u,v)])))

    return main_edge_type


def getNodeAnnotation(n, pathswithnode, pathway_colors=None, chemical=None, attr_desc=None):
    '''
    Converts the node data html for the node popup.

    :param data: dictionary of data from the NetworkX node.
    :systematic_name: The systematic name of the yeast gene 

    :returns: HTML string.
    '''

    htmlstring = ''
    uniproturl = 'http://www.uniprot.org/uniprot/%s' 
    entrezurl = 'http://www.ncbi.nlm.nih.gov/gene/%s'

    # if this is a family node:
    if len(n.split(',')) > 1:
        uids = n.split(',')
        gene_names = uniprot_to_gene[n].split(',')
        htmlstring += '<b>Gene names</b>: %s <br>' % (','.join(sorted(gene_names)))
        htmlstring += '<b>Uniprot IDs</b>: %s <br>' % (n)
        htmlstring += "<b>Uniprot IDs and Protein Names</b>:<br> <ul>"
        for gene in sorted(gene_names):
            uid = uids[gene_names.index(gene)]
            # make a bulleted list where each item in the list is:  gene name: uniprotid_link Full Protein Name
            # TODO add a link to the gene name as well
            e = db.map_id(uid, 'uniprotkb', 'GeneID')
            if len(e) > 0:
                e = e.pop() 
                entrezlink = '<a style="color:blue" href="%s" target="EntrezGene">%s</a>' % (entrezurl%e, gene) 
            else:
                entrezlink = gene
            protein_name = db.get_description(uid)
            uniprot_link = '<a style="color:blue" href="%s" target="UniProtKB">%s</a>' % (uniproturl%uid, uid)
            htmlstring += "<li>%s: %s %s </li>" % (entrezlink, uniprot_link, protein_name)
        htmlstring+='</ul>'
    else:
        #List Uniprot accession number
        uid = n
        htmlstring += '<b>Uniprot ID</b>: <a style="color:blue" href="%s" target="UniProtKB">%s</a><br>' % (uniproturl%uid, uid)
        htmlstring += '<li><i>Recommended name</i>: %s</li>' % (db.get_description(uid))
        # get the alternate names from the AltName-Full namespace
        alt_names = db.map_id(uid, 'uniprotkb', 'AltName-Full')
        if len(alt_names) > 0:
            htmlstring += '<li><i>Alternate name(s)</i>: <ul><li>%s</li></ul></li>' %('</li><li>'.join(sorted(alt_names)))

        #List EntrezGene IDs:
        for e in db.map_id(uid, 'uniprotkb', 'GeneID'):
            htmlstring += '<br><b>Entrez ID:</b> <a style="color:blue" href="%s" target="EntrezGene">%s</a>' % (entrezurl%e, e)
#
#    #List IHOP link
#    ihopurl = 'http://www.ihop-net.org/UniPub/iHOP/?search=%s&field=UNIPROT__AC&ncbi_tax_id=9606' % (uid)
#    htmlstring += '<br><b>Search in </b> <a style="color:blue" href="%s" target="IHOP">IHOP</a>' % (ihopurl)

    htmlstring += '<hr />'

    htmlstring += '<b>Paths</b>: %s<br>' %(','.join(str(i) for i in sorted(pathswithnode)))

    if attr_desc is not None and n in attr_desc:
        # now add all of the specified node annotatinos
        for attr, desc in attr_desc[n].items():
            htmlstring += "<b>%s</b>: %s</br>" % (attr, desc)

    # if the pathways are specified for this node, add them to the list
    if pathway_colors is not None:
        htmlstring += "<hr /><b>Pathways</b>:<ul>"
        for pathway in pathway_colors:
            if uid in pathway_colors[pathway]['prots']:
                pathway_link = '<a style="color:%s" href="%s">%s</a>' % (pathway_colors[pathway]['color'], pathway_colors[pathway]['link'], pathway)
                htmlstring += '<li>%s</li>' % (pathway_link)
        htmlstring+='</ul>'

    # TODO add extra node annotations
    # Chemical specific info
    if chemical is not None and chemical in CHEMICALS.chemicals and \
            n in CHEMICALS.chemicals[chemical]['accs']:
        htmlstring += "<hr /><b>Assays</b>:<ul>"
        for assay in CHEMICALS.chemicals[chemical]['accs'][n]:
            htmlstring += '<li>%s<ul>'%assay
            htmlstring += '<li>hit: <i>%d</i></li>'%CHEMICALS.chemicals[chemical]['accs'][n][assay]['hit']
            htmlstring += '<li>zscore: <i>%s</i></li>'%CHEMICALS.chemicals[chemical]['accs'][n][assay]['zscore']
            htmlstring += '<li>intended_target_type_sub: <i>%s</i></li>'%CHEMICALS.chemicals[chemical]['accs'][n][assay]['intended_target_type_sub']
            htmlstring += '<li>intended_target_family: <i>%s</i></li>'%CHEMICALS.chemicals[chemical]['accs'][n][assay]['intended_target_family']
            htmlstring += '</ul></li>'
        htmlstring+='</ul>'
    return htmlstring


##########################################################
def getEdgeAnnotation(t, h, k, evidence, family_ppi_evidence=None, attr_desc=None):
    #if t == "Q13315" and h == "Q01094":
    #    pdb.set_trace()
    annotation = '' 
    annotation+='<b>%s - %s</b></br>'%(','.join(sorted(uniprot_to_gene[t].split(','))), ','.join(sorted(uniprot_to_gene[h].split(','))))
    annotation+='<b>%s - %s</b></br>'%(t,h)
    annotation+='<b>Weight</b>: %.3f</br>' % (PPIWEIGHTS[(t,h)])
    annotation+='<b>Edge Ranking</b>: %s</br>' % (k)

    edge_str = "%s-%s" % (t,h)
    if attr_desc is not None and edge_str in attr_desc:
        # now add all of the specified edge annotatinos
        for attr, desc in attr_desc[edge_str].items():
            annotation += "<b>%s</b>: %s</br>" % (attr, desc)

    family_edge = True if len(t.split(',')) > 1 or len(h.split(',')) > 1 else False

    if family_ppi_evidence is not None and family_edge is True:
        annotation += '<hr /><h><b>Direct Sources of Evidence</b></h>'
        annotation += evidenceToHTML(t,h,evidence[(t,h)])
        if (t,h) in family_ppi_evidence:
            annotation += '<hr /><h><b>Sources of Evidence</b></h>'
            annotation += evidenceToHTML(t,h,family_ppi_evidence[(t,h)])
    else:
        annotation += '<hr /><h><b>Sources of Evidence</b></h>'
        annotation += evidenceToHTML(t,h,evidence[(t,h)])

    return annotation


##########################################################
def evidenceToHTML(t,h,evidencelist):
    annotation = '<dl>'
    sources = sorted(evidencelist)
    # move Reactome to the end as their interactions are most often ignored
    if 'Reactome' in sources:
        del sources[sources.index('Reactome')]
        sources.append('Reactome')
    if 'Reactome-FIs' in sources:
        del sources[sources.index('Reactome-FIs')]
        sources.append('Reactome-FIs')

    for source in sources:
        annotation += '<dt>%s</dt>' % (source)
        # TODO add interaction type color
        for interactiontype in evidencelist[source]:
            if interactiontype != '' and interactiontype != "None":
                # use bull instead of li to save on white space
                annotation += '&bull;&nbsp&nbsp%s <br>' % interactiontype
            for detectionmethod in evidencelist[source][interactiontype]:
                annotation += '&nbsp&nbsp&nbsp&nbsp&bull;&nbsp&nbsp'
                if detectionmethod != '':
                    filtered = False
                    for m in METHODS_TO_IGNORE:  # CSBDB
                        if m in detectionmethod:
                            filtered = True
                            break
                    # add detection method filter (i.e. gray out detection methods that we're ignoring)
                    if filtered:
                        annotation += '<FONT COLOR="gray">%s</FONT>  ' % (detectionmethod)
                    else:
                        annotation += '%s  ' % detectionmethod

                # now add the pubmed IDs. &nbsp is the html for a non-breaking space
                pub_ids = evidencelist[source][interactiontype][detectionmethod]
                #KEGG doesn't have pub ids. It has a pathway map and entry (evidence)
                try:
                    # How do we want to sort the pmid, imex, doi and MINT and such? pmid first?
                    pubmed_ids = [pubmed_id for pubmed_id in pub_ids if pubmed_id.split(':')[0] == 'pubmed' and 'None' not in pubmed_id]
                    # just sort the pubmed_ids and put them first
                    pubmed_ids = sortPubs(pubmed_ids)
                    # add the rest of the ids
                    pub_ids = pubmed_ids + [other_id for other_id in pub_ids if other_id.split(':')[0] != 'pubmed']
                    # now get the html for each of the links
                    pub_ids = [parseCSBDBpubs(pub_id) for pub_id in pub_ids if parseCSBDBpubs(pub_id) != '']
                except ValueError:
                    print("ERROR when parsing pub_ids from:")
                    print(t,h, source, interactiontype, detectionmethod, pub_ids)
                    raise

                # use a non-breaking space with a comma so they all stay on the same line
                annotation += ',&nbsp'.join(pub_ids)
                annotation += "<br>"
        annotation += '</li></ul><br>'

    return annotation


##########################################################
def parseCSBDBpubs(publication_id):
    row = publication_id.split(':')
    if row[0] == 'pubmed':
        pubmedurl = 'http://www.ncbi.nlm.nih.gov/pubmed/%s' % (row[1])
        desc = '<a style="color:blue" href="%s" target="PubMed">pmid:%s</a>' % (pubmedurl,row[1])
    elif row[0] == 'doi':
        doiurl = 'http://dx.doi.org/%s' % (row[1])
        # replace the dash with a non-linebreaking hyphen
        desc = '<a style="color:blue" href="%s" target="DOI">doi:%s</a>' % (doiurl,row[1].replace('-','&#8209;'))
    elif row[0] == 'omim':
        omimurl = 'http://omim.org/entry/%d' % (int(row[1]))
        desc = '<a style="color:blue" href="%s" target="OMIM">omim:%s</a>' % (omimurl,row[1])
    elif row[0] == 'imex':
        imexurl = 'http://www.ebi.ac.uk/intact/imex/main.xhtml?query=%s' % (row[1])
        # replace the dash with a non-linebreaking hyphen
        desc = '<a style="color:blue" href="%s" target="IMEX">imex:%s</a>' % (imexurl,row[1].replace('-','&#8209;'))
    elif row[0] == 'phosphosite':
        phosphourl = 'http://www.phosphosite.org/siteAction.action?id=%s' % (row[1])
        desc = '<a style="color:blue" href="%s" target="PSP">phosphosite:%s</a>' % (phosphourl,row[1])
    elif row[0] == 'mint':
        # MINT links are taking too long to load. The website often times out.
        # MINT links are often accompanied by IntAct links which I think have the same info
        # http://mint.bio.uniroma2.it/mint/search/interaction.do?ac=MINT-8200651
        desc = ''
    elif row[0] == 'kegg':
        # links to KEGG pathway map
        kegg_map_link = 'http://www.kegg.jp/kegg-bin/show_pathway?'
        # links to KEGG pathway entry (evidence)
        kegg_entry_link = 'http://www.kegg.jp/dbget-bin/www_bget?pathway+'
        pathway_map_link = '<a style="color:blue" href="%s%s" target="KEGG">map</a>' % (kegg_map_link,row[1])
        pathway_entry_link = '<a style="color:blue" href="%s%s" target="KEGG">evidence</a>' % (kegg_entry_link,row[1])
        desc = "%s,&nbsp%s"%(pathway_map_link,pathway_entry_link)
    elif row[0] == 'netpath':
        # TODO we can add a link to the netpath page
        desc = ''
    else:
        print("%s publication id does not have a link" % publication_id)
        # replace the dash with a non-linebreaking hyphen
        desc = publication_id.replace('-','&#8209;')

    return desc


# sort pubmed IDs based on pubmed:2334
def sortPubs(pub_ids):
    try:
        pub_ids = sorted([(pub_id.split(':')[0], int(pub_id.split(':')[1])) for pub_id in pub_ids])
        pub_ids = [':'.join([pub_id[0], str(pub_id[1])]) for pub_id in pub_ids]
    except:
        print("Warning: Could not find pubmed:",pub_ids)
    return pub_ids


##########################################################
def getGraphDescription(predfile, ppi, pathway_colors=None, chemicalID=''):
    desc = ''
    if chemicalID:
        desc += '<b>Chemical DSS ID</b>: %s <br>' % (chemicalID)
    # I used this website to get the shapes: http://shapecatcher.com/index.html
    unicode_shapes = {'triangle': '&#x25b2', 'square': '&#x25a0', 'circle': '&#x25cf',
            'arrow': '&#10230', 'T': '&#x22a3', 'dash': '&#x2014'}
    # node description table
    desc += "<b>Legend:</b>"
    desc += "<table border=\"1\">" + \
            "<tr><th>Node Style</th><th>Node Description</th></tr>" + \
            "<tr><td  align=\"center\"><font color=\"%s\" size=\"8\">%s;</font></td><td>Source Receptor</td></tr>" % (NODE_COLORS['source'], unicode_shapes['triangle']) + \
            "<tr><td  align=\"center\"><font color=\"%s\" size=\"8\">%s;</font></td><td>Target TF</td></tr>" % (NODE_COLORS['target'], unicode_shapes['square']) + \
            "<tr><td  align=\"center\"><font color=\"%s\" size=\"8\">%s;</font></td><td>Receptor</td></tr>" % (NODE_COLORS['intermediate_rec'], unicode_shapes['triangle']) + \
            "<tr><td  align=\"center\"><font color=\"%s\" size=\"8\">%s;</font></td><td>TF</td></tr>" % (NODE_COLORS['intermediate_tf'], unicode_shapes['square']) + \
            "<tr><td  align=\"center\"><font color=\"%s\" size=\"8\">%s;</font></td><td>Intermediate Protein</td></tr>" % (NODE_COLORS['default'], unicode_shapes['circle']) + \
            "</table>"
    desc += "<br><b>Source Receptor</b>: Receptor used as source for KSP"
    desc += "<br><b>Target TF</b>: Transcription Factor/Regulator used as target for KSP"
    # also add the edge style and description table
    # 'physical' 'phosphorylation' 'enzymatic' 'activation' 'inhibition' 'spike_regulation'
 
    desc += "<br><br><table border=\"1\">" + \
            "<tr><th>Edge Style</th><th>Edge Description</th></tr>" + \
            "<tr><td>&nbsp;&nbsp;<b><font color=\"%s\" size=\"5\">%s;</font></b></td>       <td>Phosphorylation</td></tr>"    % (EDGE_COLORS['phosphorylation'], unicode_shapes['arrow']) + \
            "<tr><td>&nbsp;&nbsp;<b><font color=\"%s\" size=\"5\">%s;</font></b></td>       <td>Enzymatic Reaction</td></tr>" % (EDGE_COLORS['enzymatic'], unicode_shapes['arrow']) + \
            "<tr><td>&nbsp;&nbsp;<b><font color=\"%s\" size=\"5\">%s;</font></b></td>       <td>KEGG Activation</td></tr>"    % (EDGE_COLORS['activation'], unicode_shapes['arrow']) + \
            "<tr><td>&nbsp;&nbsp;&nbsp;<b><font color=\"%s\" size=\"5\">%s;</font></b></td> <td>KEGG Inhibition</td></tr>"    % (EDGE_COLORS['inhibition'], unicode_shapes['T']) + \
            "<tr><td>&nbsp;&nbsp;<b><font color=\"%s\" size=\"5\">%s;</font></b></td>       <td>SPIKE regulation</td></tr>"   % (EDGE_COLORS['spike_regulation'], unicode_shapes['arrow']) + \
            "<tr><td>&nbsp;&nbsp;&nbsp;<b><font color=\"%s\" size=\"5\">%s;</font></b></td> <td>Physical</td></tr>"           % (EDGE_COLORS['physical'], unicode_shapes['dash']) + \
            "</table>"

    #desc += 'Edges can have multiple sources of evidence and therefore multiple interaction types. Edge style was given based on the order in the table.'
    desc += '<br>Enzymatic reactions are reactions such as %s, etc.' % (", ".join(sorted(["dephosphorylation", "ubiquitination", "methylation", "glycosylation"])))
    desc += '<br>Phosphorylation is separate because it occurs more frequently.'

    desc += '<hr /><b>Sources of Evidence</b><br> Sources are black if they were used to construct the interactome ' + \
            'and <FONT COLOR="gray">gray</FONT> if they were not used to construct the interactome.' 
    desc +='<ul>'
    # main databases used to construct the interactome
    for db in ['BioGrid', 'DIP', 'KEGG', 'IntAct', 'Mint', 'NetPath', 'Phosphosite Plus', 'Spike']:
        desc+='<li><FONT COLOR="black">%s</FONT></li>' % db
    # Reactome may have some interactions that weren't filtered, but most are
    for db in ['Reactome', 'Reactome-FIs']:
        desc+='<li><FONT COLOR="grey">%s</FONT></li>' % db
    desc+='</ul>'

    desc += '<hr /><b>Path/Edges file</b>: %s<br>' % (predfile)
    desc += '<hr /><b>PPI background interactome file</b>: %s<br>' % (ppi)

    # TODO add descriptions of the pathway colors
    if pathway_colors is not None:
        #print("Adding color description from %s " % (pathway_colors))
        table =  "<br><br><table border=\"1\">"
        table += "<tr><th>Pathway</th><th>Estimate</th></tr>"
        for pathway in pathway_colors:
            pathway_link = '<a style="color:black" href="%s" target="KEGG">%s</a>' % (pathway_colors[pathway]['link'], pathway)
            table += "<tr><td bgcolor=\"%s\">%s</td> <td>%s</td></tr>" % (pathway_colors[pathway]['color'], pathway_link, pathway_colors[pathway]['estimate'])
        table += "</table>"
        desc += table


    return desc


def parseArgs(args):
    ## Parse command line args.
    usage = 'post_to_graphspace_human.py [options]\n'
    parser = OptionParser(usage=usage)
    parser.add_option('', '--paths', type='string', metavar='STR',
                      help='file of Pathlinker path results. Use this or the --ranked-edges option')
    parser.add_option('', '--ranked-edges', type='string', metavar='STR',
                      help='file of Pathlinker ranked edges results. Use this or the --paths option')
    parser.add_option('', '--ppi', type='string', metavar='STR',
                      help='file of weighted directed edges. Tab-delimited file has columns TAIL,  HEAD,  WEIGHT,  EVIDENCE.  WEIGHT must be nonnegative. Required')
    parser.add_option('', '--mapping-file', type='string', metavar='STR',
                      help='File used to map to a different namespace. Network/edge IDs (uniprot ids) should be in the first column with the other namespace (gene name) in the second')
    parser.add_option('', '--datadir', type='str', metavar='STR', default='/home/jeffl/svnrepo/data',
                      help='Data directory from SVN. Required.')
    parser.add_option('','--version',type='string',default='2016_05pathlinker',
                      help='Version of the evidence file, which contains all the compiled evidence for all of the interactions. ')  #Also used to get the uniprot_to_gene mapping. Options are %s.  Default is "2016_05pathlinker".' % (','.join(getev.VERSIONS)))
    parser.add_option('','--ev-file',type='string',
                      help='Path to evidence file if a different file from the default is to be used.')
    #parser.add_option('', '--chemicalID', type='string', metavar='STR',
    #                  help='Chemical ID of interest. Only the actual ID (i.e. DSSTox_GSID_47628 = 47628).' +
    #                  'n\tIf this option is specified,  the non-hit proteins will have double black border and hit proteins have a double brown border')
    parser.add_option('', '--sourcetarget', type='string', metavar='STR',
                      help='Source-target file used as input to pathlinker. Required.')
    parser.add_option('', '--human-rec-tfs', type='string', metavar='STR',
                      help='List of human receptors and tfs. Useful to add shape (triangle or square) to intermediate nodes. Default=\'/data/jeff-law/projects/2016-02-toxcast/inputs/pathlinker-data/human-rec-tfs.txt\'')

    # posting options
    parser.add_option('-U', '--username', type='string', metavar='STR', 
                      help='GraphSpace account username to post graph to. Required')
    parser.add_option('-P', '--password', type='string', metavar='STR', 
                      help='Username\'s GraphSpace account password. Required')
    parser.add_option('', '--graph-name', type='string', metavar='STR', default='test',
                      help='Graph name for posting to GraphSpace. Default = "test".')
    parser.add_option('', '--outprefix', type='string', metavar='STR', default='test',
                      help='Prefix of name to place output files. Required.')
    parser.add_option('', '--group', type='string', metavar='STR',
                      help='Group to share the graph with.')
    parser.add_option('', '--group-owner', type='string', metavar='STR', default='jeffl@vt.edu',
                      help='Owner of the group.')
    parser.add_option('', '--tag', action="append", type='string', metavar='STR',
                      #help='Tag to put on the graph.')
                      help='Tag to put on the graph. Can list multiple tags (for example --tag tag1 --tag tag2)')
    parser.add_option('-k', '--k-limit', type='int', default=1000,
                      help='The k limit of paths to post to graphspace.')
    parser.add_option('', '--apply-layout', type='string', metavar='STR', 
                      help='Specify the name of a graph from which to apply a layout. Layout name specified by the --layout-name option.')
    parser.add_option('', '--layout-name', type='string', metavar='STR', default='layout1',
                      help="Name of the layout of the graph specified by the --apply-layout option to apply. Default: 'layout1'")

    # extra options
    parser.add_option('', '--graph-attr', type='string', metavar='STR',
            help='Tab-separated File used to specify graph attributes 1st column: style, 2nd: style attribute, 3rd: list of uniprot ids (nodes/edges) separated by |, 4th: Description.')
    parser.add_option('', '--casestudy', action="store_true", default=False,
                      help="Use the Case study colors and labels (no k in labels, gray border around nodes, all directed interactions orange)")
    parser.add_option('', '--include-parent-nodes', action="store_true",
                      help="Include source, target, intermediate parent nodes, and whatever parent nodes are in the --graph-attr file")
    parser.add_option('', '--node-overlap', type='string',
                      help="File containing the # of paths each node participates in")
    parser.add_option('', '--edge-overlap', type='string',
                      help="File containing the # of paths each edge participates in")
    #parser.add_option('', '--min-num-networks', type='int',
    #                  help="Minimum number of networks an edge needs to be in to be included in the graph")
    parser.add_option('', '--ctd-support', type='string',
                      help="Color nodes according to the phosphorylation support in CTD chemical-gene-interactions file.")
    parser.add_option('', '--sig-chemicals', type='string',
                      help="File containing the list of significant chemicals.")
    parser.add_option('', '--overlap-cutoff', type='float', default=0.5,
                      help="Float representing the % of networks nodes/edges must be in to be included")

    (opts, args) = parser.parse_args(args)

    #if not opts.ranked_edges and not opts.paths:
    #    parser.print_help()
    #    sys.exit("\n--paths or --ranked_edges option required")
    #if not opts.ppifile:
    #    parser.print_help()
    #    sys.exit("\n--ppifile option required")
    if not opts.sourcetarget:
        parser.print_help()
        sys.exit("\n--sourcetarget option required")

    return opts, args


if __name__ == '__main__':
    main(sys.argv)
