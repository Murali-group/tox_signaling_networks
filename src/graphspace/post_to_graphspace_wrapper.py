
# code to post toxcast results to graphspace

import argparse
import os
import sys
import json
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from src.utils import file_utils as utils
#from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
# used for coloring the nodes according to hit/nonhit values
#sys.stderr.write("ImportError: toxcast_utils.py not found")

from src import toxcast_utils as t_utils
from src import toxcast_settings as t_settings
#import post_to_new_graphspace_evidence as post_to_gs
from src.graphspace import post_to_graphspace_base as gs_base
from src.graphspace import gs_utils

toxcast_data = t_utils.parse_toxcast_data.ToxCastData() 
chemIDtoName, chemNametoID = t_utils.getChemicalNameMaps(toxcast_data)
chemIDtoCAS, chemCAStoID = t_utils.getChemicalCASMap(toxcast_data)

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
        'No Enriched Terms': '#D8D8D8',  # gray
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
        'enzymatic': 'brown',
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

# define a set of colors to pull from
# TODO what if there are more than 16 pathways?
# here's a great set of colors: http://clrs.cc/
COLORS = [
    # light-green, lime-green?, turquoise, aqua, purple, 
    "#56d67a", "#00cc99", "#0099cc", "#7FDBFF", "#9966ff",
    # orange-red, orange, light-green, yellow, blue,
    "#ff5050", "#e67300", "#86b300", "#FFDC00", "#0074D9",
    # blue2,    pink,      orange2,  maroon, light-green2, light_purple, dark green
    "#6699ff", "#f7608b", "#f79760", "#85144b", "#4efc77", "#eb46e2", "#34a352",
    # LIGHT purple, DARK purple, dark brown, dark blue, salmon, orange
    "#cc99ff", "#6600cc", "#996633", "#006699", "#ffcccc", "#e68a00",
]

# for some of the terms that are common for multiple chemicals, give them the same color
FIXED_TERM_COLORS = {
    "GO:0042493": "#85144b",  # response to drug
    "GO:0043066": "#f79760",  # neg. reg. apoptotic process
    "GO:0016032": "#f7608b",  # viral process
    "GO:0038095": "#6699ff",  # Fc-epsilon receptor\nsignaling pathway
    "GO:0030168": "#0099cc",  # platelt activation
    # make this term red because it makes sense
    "GO:1900034": "#ff5050",  # regulation of cellular response to heat
    }


def main(**kwargs):
    """
    *opts*: the optparse options to pass to post_to_new_graphspace_evidence
    *kwargs*: the kwargs version of opts
    """
    #chemicals = kwargs.get('chemicals')
    #version = kwargs['version']
    #del kwargs['version']
    #if chemicals is None:
        #chemicals = sorted(readItemList("%s/chemicals.txt" % (t_settings.INPUTSPREFIX)))
    call_post_to_graphspace(**kwargs)


def write_revigo_color_files(chemicals, RESULTSPREFIX, **kwargs):
    """
    If a file downloaded from REVIGO is passed in, then use that to set the term colors and boxes.
    Otherwise, just remove the most frequent term. Hopefully there isn't too much overlap in the remaining terms.
        TODO Another possible strategy is to cluster the terms myself and select a single term per cluster.
    """

    # assign a color to each term
    out_dir = "%s/graphspace/colors" % (RESULTSPREFIX)
    t_utils.checkDir(out_dir)
    #print("Writing REVIGO colors to %s for %d chemicals. (limit of %d colors)" % (out_dir, len(chemicals), len(colors)))

    chem_color_files = {}

    for chemical in chemicals:
        out_prefix = "%s/%s" % (out_dir, chemical)
        out_file = "%s-colors.tsv" % (out_prefix) 
        # first read the david results file
        david_file = "%s/stats/go-analysis/chemical-reports/%s.txt" % (RESULTSPREFIX, chemical)
        print("reading %s" % (david_file))
        df = pd.read_csv(david_file, sep='\t')
        print(df.head())
        # get just the term ids
        df['term'] = df['Term'].apply(lambda x: x[:x.find('~')])
        df['name'] = df['Term'].apply(lambda x: x.split('~')[1])
        # build a dictionary from the term to the prots
        orig_term_names = dict(zip(df['term'], df['name']))
        term_prots = {t: prots.replace(', ', '|') for t, prots in zip(df['term'], df['Genes'])}
        term_pvals = dict(zip(df['term'], df['Bonferroni']))
        term_names = dict(zip(df['term'], df['name']))
        name_to_term = dict(zip(df['name'], df['term']))

        # read the revigo file and extract the GO term info
        if kwargs.get('revigo_file'): 
            if not os.path.isfile(kwargs['revigo_file']):
                print("ERROR: --revigo-file '%s' not found." % (kwargs['revigo_file']))
                sys.exit()
            print("reading %s" % (kwargs['revigo_file']))
            df_r = pd.read_csv(kwargs['revigo_file'], sep=',')
            print(df_r.head())
            # sort by pval
            #df_r = df_r.sort_values("log10 p-value")
            term_names = dict(zip(df_r['term_ID'], df_r['description']))
            selected_terms = list(term_names.keys())
        elif kwargs.get('term_counts_file'):
            if not os.path.isfile(kwargs['term_counts_file']):
                print("ERROR: --term-counts-file '%s' not found." % (kwargs['term_counts_file']))
                sys.exit()
            print("reading %s" % (kwargs['term_counts_file']))
            df_r = pd.read_csv(kwargs['term_counts_file'], sep='\t', names=['term_name', 'count'])
            print(df_r.head())
            freq_cutoff = kwargs.get('freq_cutoff', .75)
            print("applying a frequency cutoff of %s" % (freq_cutoff))
            df_r['freq'] = df_r['count'] / df_r['count'].max() 
            term_freq = dict(zip(df_r['term_name'], df_r['freq']))
            # sort by pval
            df = df.sort_values("Bonferroni")
            # apply a cutoff of 0.01
            df = df[df['Bonferroni'] < 0.01]
            selected_terms = []
            for name in df['name']:
                if term_freq[name] < freq_cutoff:
                    selected_terms.append(name_to_term[name])

        term_popups = {}
        link_template = "<a style=\"color:blue\" href=\"https://www.ebi.ac.uk/QuickGO/GTerm?id=%s\" target=\"DB\">%s</a>"
        for term in selected_terms:
            term_link = link_template % (term, term)
            popup = "<b>QuickGO</b>: %s" % (term_link)
            popup += "<br><b>p-value</b>: %0.2e" % (float(term_pvals[term]))
            term_popups[term] = popup

        function_colors = write_colors_file(out_file, selected_terms, term_names, term_prots, term_popups)
        chem_color_files[chemical] = out_file 

        new_func_colors = defaultdict(dict)
        for term in function_colors:
            new_func_colors[term]['prots'] = term_prots[term]
            new_func_colors[term]['color'] = function_colors[term]
            new_func_colors[term]['link'] = "https://www.ebi.ac.uk/QuickGO/GTerm?id=%s" % (term)
            new_func_colors[term]['name'] = orig_term_names[term]
            # if uid in pathway_colors[pathway]['prots']:
            #     pathway_link = '<a style="color:%s" href="%s">%s</a>' % (pathway_colors[pathway]['color'], pathway_colors[pathway]['link'], pathway)
    return chem_color_files, new_func_colors


def write_colors_file(out_file, functions, function_names, function_prots, function_popups=None, colors=None):
    # first, shorten the funciton names and make them wrap over two lines
    function_names = shorten_names(function_names)
    if colors is None:
        colors = COLORS
    if len(functions) > len(colors):
        print("\tWarning: # functions %d exceeds # colors %d. Limiting functions to the %d colors" % (len(functions), len(colors), len(colors)))
        prots = set([p for prots in function_prots.values() for p in prots.split('|')])
        prots_used = set()
        new_funcs = []
        #for f in sorted(function_prots, key=lambda f: len(function_prots[f]), reverse=True):
        for f in functions:
            f_prots = set(function_prots[f].split('|'))
            # only keep functions if they have some unique annotations
            if len(f_prots - prots_used) < 2:
                continue
            else:
                new_funcs.append(f)
                prots_used.update(f_prots)
        functions = new_funcs
        if len(new_funcs) > len(colors):
            # limit the pathways to the number of colors
            functions = new_funcs[:len(colors)]

    function_colors = {}
    available_colors = set(colors)
    for function in functions:
        if function in FIXED_TERM_COLORS:
            function_colors[function] = FIXED_TERM_COLORS[function]
            available_colors.discard(function_colors[function])
    available_colors = [colors[i] for i, color in enumerate(colors) if color in available_colors]
    for i, function in enumerate(sorted(set(functions) - set(function_colors.keys()))):
        function_colors[function] = available_colors[i] 

    # if no popup was provided, just show an empty popup
    if function_popups is None:
        function_popups = {function:"" for function in functions} 

    print("\tWriting graphspace colors file: %s" % (out_file))
    ## UPDATE 2017-07-13: create a compound or "parent" node for each goterm where the name of the parent node is the goterm name
    ## TODO The 4th column is the description/popup of the parent node
    # write the prots along with a color for posting to graphspace
    with open(out_file, 'w') as out:
        # first write the style of the individual nodes
        #out.write("#style\tstyle_val\tprots\tdescription\n")
        out.write('\n'.join(['\t'.join(['color', function_colors[function], function_prots[function], '-']) for function in functions])+'\n')
        # then write the parent nodes
        #out.write("#parent\tfunction\tprots\tpopup\n")
        for function in functions:
            out.write('\t'.join(['parent', function_names[function], function_prots[function], function_popups[function]])+'\n')
        # then write the style of the parent node
        #out.write("#color\tcolor_val\tparent_node\tdescription\n")
        out.write('\n'.join(['\t'.join(['color', function_colors[function], function_names[function], '-']) for function in functions])+'\n')
    return function_colors


def shorten_names(function_names):
    new_func_names = {}
    for term, name in function_names.items():
        if len(name) > 20:
            split_name = name.split(' ')
            curr_len = 0
            curr_name = ''
            for i, word in enumerate(split_name):
                curr_name += word
                curr_len += len(word) + 1
                if curr_len > 20 and i < len(split_name)-1:
                    curr_name += '\\n'
                    curr_len = 0 
                else:
                    curr_name += ' '
                if len(curr_name) > 60:
                    break
            name = curr_name
        new_func_names[term] = name
    return new_func_names


def get_ctd_support(chemical, prednodes, ctd_support_file):
    print("Getting CTD support counts from %s" % (ctd_support_file))
    num_intxs_per_gene = {}
    num_intxs_per_node = {}
    for cas, gene, interaction_action, pubmedids in utils.readColumns(ctd_support_file, 3, 4, 10, 11):
        if chemIDtoCAS[chemical] != cas:
            continue
        if 'phosphorylation' in interaction_action:
            if gene not in num_intxs_per_gene:
                num_intxs_per_gene[gene] = 0
            num_intxs_per_gene[gene] += 1 

    for n in prednodes:
        gene = uniprot_to_gene[n]
        if gene in num_intxs_per_gene:
            # for now, take the node in the family node with the maximum support
            num_intxs_per_node[n] = num_intxs_per_gene[gene]
    print("\nOf the %d prots with (de)phosphorylation evidence in CTD, %d of the %d net nodes overlap" % (len(num_intxs_per_gene), len(num_intxs_per_node), len(prednodes)))

    # # write these counts to a table
    # out_file = "%s-ctd-support.txt" % (opts.node_overlap)
    # print("Writing support counts to %s" % (out_file))
    # with open(out_file, 'w') as out:
    #     out.write("#uniprot\tgene\tmax_num_phospho_intxs\n")
    #     # write the sorted results to a file
    #     out.write('\n'.join(["%s\t%s\t%d" % (N, uniprot_to_gene[N], num_intxs_per_node[N]) for N in sorted(num_intxs_per_node, key=num_intxs_per_node.get, reverse=True)]) + '\n')
    return num_intxs_per_node


def call_post_to_graphspace(version, chemicals, **kwargs):
    INPUTSPREFIX, RESULTSPREFIX, interactome = t_settings.set_version(version)
    t_utils.checkDir("%s/graphspace" % (RESULTSPREFIX))
    # write the color files of each chemical and return a dictionary of the chemical and its color file
    #if kwargs['revigo_colors']:
    chemical_color_files = None 
    if kwargs.get('revigo_file') or kwargs.get('term_counts_file'):
        chemical_color_files, function_colors = write_revigo_color_files(chemicals, RESULTSPREFIX, forced=kwargs['forcepost'], **kwargs)
        kwargs['function_colors'] = function_colors
        print(chemical_color_files)
    kwargs['tags'] = kwargs['tags'] + [version] if kwargs.get('tags') else [version]
    # post everything to graphspace!
    for chemical in chemicals:
        # get the chemical name. make sure it doesn't have any '%' or '[',']' as that will break posting
        chemName = chemIDtoName[chemical].replace('%','').replace('[','').replace(']','')
        rec_tfs_file = t_settings.REC_TFS_FILE % (INPUTSPREFIX, chemical)
        cyclinker_output_file = "%s/cyclinker/%s-paths.txt"%(RESULTSPREFIX, chemical)
        output_json = "%s/graphspace/%s-graph%s.json" % (RESULTSPREFIX, chemical, kwargs.get('name_postfix',''))
        proteins, num_paths = t_utils.getProteins(paths=cyclinker_output_file, max_k=kwargs['k_to_post'], ties=True)
        if not kwargs['forcepost'] and os.path.isfile(output_json): 
            print("%s already exists. Use --forcepost to overwrite it" % (output_json))
        else:
            # TEMP: The evidence file I used is different from the default:
            ev_file = "/data/jeff-law/svnrepo/data/interactions/compiled/2018_01/with-d2d/2018_01pathlinker-nokegg.tsv"
            build_graph_and_post(
                version, interactome, rec_tfs_file, RESULTSPREFIX, cyclinker_output_file,
                chemical, max_k=num_paths, graph_name="%s-%s-%s%s"%(chemName,chemical,version,kwargs.get('name_postfix','')),
                #name_postfix='-'+version, tag=version, chemical_color_file=)
                graph_attr_file=chemical_color_files.get(chemical) if chemical_color_files is not None else None,
                ev_file=ev_file, out_pref="%s/graphspace/%s%s" % (RESULTSPREFIX, chemical,kwargs.get('name_postfix','')),
                **kwargs)


def build_graph_and_post(
        version, interactome, rec_tfs_file, RESULTSPREFIX, paths_file,
        chemical, max_k=200, graph_name="test",
        #postfix='-'+version, tag=version, chemical_color_file=)
        graph_attr_file=None, ev_file=None,
        datadir="/home/jeffl/svnrepo/data", **kwargs):
    # get the "evidence version" which is used to get the CSBDB related files
    ev_version = t_utils.get_ev_version(version)

    PPI = interactome
    lines = utils.readColumns(PPI,1,2,3)
    global PPIEDGES, PPIWEIGHTS
    PPIEDGES = [(u,v) for u,v,w in lines]
    PPIWEIGHTS = {(u,v):float(w) for u,v,w in lines}

    prededges = readNetwork(paths=paths_file, k_limit=max_k)

    sources = set()
    targets = set()
    lines = utils.readColumns(rec_tfs_file,1,2)
    sources = set([acc for acc, node_type in lines if node_type.lower() in ['source', 'receptor']])
    targets = set([acc for acc, node_type in lines if node_type.lower() in ['target', 'tf']])
    # human_rec = set() 
    # human_tfs = set()
    # if opts.human_rec_tfs is not None:
    #     # also get the human rec and tfs
    #     lines = utils.readColumns(opts.human_rec_tfs,1,2)
    #     human_rec = set([acc for acc, node_type in lines if node_type.lower() in ['source', 'receptor']])
    #     human_tfs = set([acc for acc, node_type in lines if node_type.lower() in ['target', 'tf']])

    # TODO tempororay fix for family nodes
    prednodes = set([t for t,h in prededges]).union(set([h for t,h in prededges]))

    global uniprot_to_gene
    #uniprot_to_gene = utils.readDict(getev.getMappingFile(ev_version, datadir), 1, 2)
    uniprot_to_gene = utils.readDict(kwargs['mapping_file'], 1, 2)

    # get attributes of nodes and edges from the graph_attr file
    graph_attr = {}
    # description of a style, style_attr tuple 
    attr_desc = {}
    if graph_attr_file is not None:
        graph_attr, attr_desc = gs_base.readGraphAttr(graph_attr_file) 

    # if opts.chemicalID:
    #     global CHEMICALS
    #     CHEMICALS = t_utils.loadChemicalMap(PPI)

    if kwargs.get('ctd_support_file'):
        num_intxs_per_node = get_ctd_support(chemical, prednodes, kwargs['ctd_support_file'])
        # set the double border attribute for nodes nodes with any CTD support
        for n in num_intxs_per_node:
            graph_attr[n]['style'] = 'double'
            graph_attr[n]['border_color'] = 'maroon'
            graph_attr[n]['border_width'] = 10

    # set the case study colors - all nodes are gray by default
    if kwargs.get('case_study') is True:
        if graph_attr_file is not None:
            # if there are other colors present, then make the sources and targets gray by default because they could have other colors
            NODE_COLORS.update(CASESTUDY_NODE_COLORS)
            EDGE_COLORS.update(CASESTUDY_EDGE_COLORS)

    # get the evidence supporting each edge
    evidence, edge_types, edge_dir = gs_utils.getEvidence(prededges.keys(), evidence_file=ev_file)

    # Now post to graphspace!
    #G = gs.constructGraph(pred_edges, node_labels=uniprot_to_gene, graph_attr=graph_attr, popups=popups)
    #G = gs_base.constructGraph(prededges, node_labels=uniprot_to_gene, graph_attr=graph_attr, attr_desc=attr_desc)
    G = constructGraph(
        prededges, sources, targets,
        node_labels=uniprot_to_gene,
        evidence=evidence, edge_types=edge_types, edge_dir=edge_dir,
        graph_attr=graph_attr, attr_desc=attr_desc, **kwargs)
    print("Graph has %d nodes and %d edges" % (G.number_of_nodes(), G.number_of_edges()))

    # put the parent nodes and the nodes in the parent nodes in a grid layout automatically
    print("Setting the x and y coordinates of each node in a grid layout")
    # relabel the nodes to their names
    graph_attr = {uniprot_to_gene.get(n,n): attr for n, attr in graph_attr.items()}
    layout = gs_utils.grid_layout(G, graph_attr)
    for node, (x,y) in layout.items():
        G.set_node_position(node_name=node, x=x, y=y)

    # before posting, see if we want to write the Graph's JSON to a file
    if kwargs.get('out_pref') is not None:
        print("Writing graph and style JSON files to:\n\t%s-graph.json \n\t%s-style.json" % (kwargs['out_pref'], kwargs['out_pref']))
        with open(kwargs['out_pref']+"-graph.json", 'w') as out:
            json.dump(G.get_graph_json(), out, indent=2)
        with open(kwargs['out_pref']+"-style.json", 'w') as out:
            json.dump(G.get_style_json(), out, indent=2)

    G.set_tags(kwargs.get('tags',[]))
    G.set_name(graph_name)

    gs_base.post_graph_to_graphspace(
            G, kwargs['username'], kwargs['password'], graph_name, 
            apply_layout=kwargs['apply_layout'], layout_name=kwargs['layout_name'],
            group=kwargs['group'], make_public=kwargs['make_public'])


def readNetwork(paths=None, ranked_edges=None, k_limit=200, no_k=False):
    """ Read the PathLinker paths or ranked_edges output. 
        Get all of the edges that have a k less than the k_limit.
    """
    if no_k is False:
        if paths is not None:
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

        if ranked_edges is not None:
            # Predicted edges from pathlinker
            lines = utils.readColumns(ranked_edges,1,2,3)
            # get all of the edges that have a k less than the k_limit
            prededges = {(u,v):int(k) for u,v,k in lines if int(k) <= k_limit}
    else:
        if ranked_edges:
            # Set of edges from another source such as a pathway
            lines = utils.readColumns(ranked_edges,1,2)
            # keep the edges in a dictionary to work with the rest of the code
            prededges = {(u,v):None for u,v in lines}

    return prededges


def constructGraph(
        prededges, sources, targets, node_labels={},
        evidence={}, edge_types={}, edge_dir={},
        graph_attr={}, attr_desc={}, **kwargs):
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
    # edges_file = kwargs['ranked_edges'] if kwargs.get('ranked_edges') is not None else kwargs['paths']
    # # get metadata
    # desc = getGraphDescription(edges_file, kwargs['ppi'], pathway_colors=None, chemicalID=kwargs['chemicalID'])
    # metadata = {'description':desc,'tags':[], 'title':''}
    # if kwargs['tag']:
    #     metadata['tags'] = [kwargs['tag']]
    # if kwargs['chemicalID']:
    #     metadata['title'] = CHEMICALS.chemDSStoName[kwargs['chemicalID']]

    # # get the evidence from get_interaction_evidence.py
    # evidence_file = kwargs.get('evidence_file') 
    # if evidence_file is None:
    #     evidence_file = getev.getEvidenceFile(kwargs['version'], kwargs['datadir'])
    # # TODO make it an option to add evidence to family edges
    # add_ev_to_family_edges = True 
    # family_ppi_evidence = None 
    # if add_ev_to_family_edges:
    #     evidence, family_ppi_evidence, edge_types, edge_dir = getev.getEvidence(prededges.keys(), evidence_file=evidence_file, split_family_nodes=True, add_ev_to_family_edges=add_ev_to_family_edges)
    # else:
    #     evidence, edge_types, edge_dir = getev.getEvidence(prededges.keys(), evidence_file=evidence_file, add_ev_to_family_edges=add_ev_to_family_edges)

    prednodes = set([t for t,h in prededges]).union(set([h for t,h in prededges]))

    # GSGraph does not allow adding multiple nodes with the same name.
    # Use this set to make sure each node has a different gene name.
    # if the gene name is the same, then use the gene + uniprot ID instead
    genes_added = set()
    # set of parent nodes to add to the graph
    parents_to_add = {}

    ## add GraphSpace/Cytoscape.js attributes to all nodes.
    for n in prednodes:
        #default is gray circle
        node_type = 'default'
        # default parent
        #parent = 'Intermediate Proteins'
        parent = 'No Enriched Terms'
        if n in sources:
            # if n is the source, make it a blue triangle
            node_type = 'source'
            # set the parent node for this source 
            parent = 'Sources'
        elif n in targets:
            # if n is a taret, make it a yellow square
            node_type = 'target'
            parent = 'Targets'
        # # also check if the individual nodes in the family node are counted as rec or TFs
        # elif n in extra_rec or len([p for p in n.split(',') if p in extra_rec]) > 0:
        #     # gray triangle
        #     node_type = 'intermediate_rec'
        #     parent = 'Intermediate Receptors'
        # # also check if the individual nodes in the family node are counted as rec or TFs
        # elif n in extra_tfs or len([p for p in n.split(',') if p in extra_tfs]) > 0:
        #     # gray square
        #     node_type = 'intermediate_tf'
        #     parent = 'Intermediate TFs'

        if n not in graph_attr:
            graph_attr[n] = {}
            graph_attr[n]['parent'] = parent
        # only overwrite the parent if it is a source or target
        elif node_type == 'source' or node_type == 'target':
            graph_attr[n]['parent'] = parent

        # set the name of the node to be the gene name and add the k to the label
        gene_name = node_labels[n]

        # if this gene name node was already added, then add another node with the name: gene-uniprot
        # TODO some uniprot IDs map to the same Gene Name (such as GNAS)
        if gene_name in genes_added:
            gene_name = "%s-%s" % (gene_name, n)
            node_labels[n] = gene_name
            #continue
        genes_added.add(gene_name)

        short_name = gene_name
        # # TODO if the the family name is too long, then just choose one of the genes and add -family to it
        # if len(gene_name) > 10:
        #     short_name = "%s-family" % (gene_name.split(',')[0])

        if kwargs.get('no_k') is True:
            pathswithnode = None 
            k_value = None
        else:
            edgeswithnode = set([(t,h) for t,h in prededges if t==n or h==n])
            pathswithnode = set([int(prededges[e]) for e in edgeswithnode])
            k_value = min(pathswithnode)
        if kwargs.get('case_study') is True:
            # set all of the sources and targets to 1 so they will all always be visible
            k_value = 0 if node_type in ["source", "target"] else k_value

        attr_dict = {}
        if kwargs.get('parent_nodes') is True:
            # set the parent if specified
            if n in graph_attr and 'parent' in graph_attr[n]:
                # set the parent of this node
                parent = graph_attr[n]['parent']
                attr_dict['parent'] = parent
            # also add this parent to the set of parent/compound nodes to add to the graph
            # keep track of the lowest k value so it will work correctly with the sliding bar
            if parent not in parents_to_add or k_value < parents_to_add[parent]:
                parents_to_add[parent] = k_value

        node_popup = buildNodePopup(n, pathswithnode, pathway_colors=kwargs.get('function_colors'), chemical=kwargs.get('chemicalID'))

        label = short_name
        if kwargs.get('case_study', False) is False and kwargs.get('no_k',False) is False:
            label="%s\n%d"%(short_name,k_value)

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
        border_color = None
        if kwargs.get('case_study') is True:
            border_color = "#7f8184"
        bubble = None
        if n in graph_attr:
            # TODO allow for any attribute to be set. Need to update add_node_style first so that attributes aren't overwritten
            #for style in graph_attr[n]:
            #    attr_dict[style] = graph_attr[n][style]
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
        border_color = color if border_color is None else border_color
        # I updated the bubble function in graphspace_python gsgraph.py so it wouldn't overwrite the border color.
        bubble = color if bubble is None else bubble
        luminance = get_color_luminance(color)
        # this sets the text color to white for dark colors
        if luminance < 0.45:
            attr_dict['color'] = 'white'

        # if kwargs['chemicalID'] and n in CHEMICALS.chemical_protein_hit[kwargs['chemicalID']] \
        #         and (kwargs.get('case_study', False) is False or (node_type == "source" or node_type == "target")):
        #     #print(n, CHEMICALS.chemicals[chemicalID]['accs'][n])
        #     style = 'double'
        #     border_width = 10
        #     if CHEMICALS.chemical_protein_hit[kwargs['chemicalID']][n] > 0:
        #         border_color = 'maroon'
        #     else: 
        #         border_color = 'black'
        # if this is a source or target, then apply the double border automatically
        #if kwargs['chemicalID'] and (node_type in ["source", "target"]): 
        # TODO need to update graphspace so the double border color applies
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
        label = label.replace("\\n","\n")

        #popup = parent
        # leave off the label for now because it is written under the edges and is difficult to see in many cases
        # I requested this feature on the cytoscapejs github repo: https://github.com/cytoscape/cytoscape.js/issues/1900
        G.add_node(parent, popup=popup, k=k_value, label=label)

        parent_label = parent.replace('\n','\\n')
        # the only style we need to set for the parent/compound nodes is the color
        # default color for sources and targets, and intermediate nodes
        if parent in NODE_COLORS:
            color = NODE_COLORS[parent]
        if parent_label in graph_attr and 'color' in graph_attr[parent_label]:
            color = graph_attr[parent_label]['color']
        attr_dict = {} 
        # set the background opacity so the background is not the same color as the nodes
        attr_dict["background-opacity"] = 0.3
        attr_dict['font-weight'] = "bolder"
        attr_dict['font-size'] = "32px"
        attr_dict['text-outline-width'] = 3
        luminance = get_color_luminance(color)
        # this sets the text color to white for dark colors
        if luminance < 0.45:
            attr_dict['color'] = 'white'
        # remove the border around the parent nodes 
        border_width = 0
        valign="bottom"

        G.add_node_style(parent, attr_dict=attr_dict, color=color, border_width=border_width, bubble=color, valign=valign)

    # Add all of the edges and their Graphspace/Cytoscape.js attributes
    for (u,v) in prededges:
        # get the main edge type
        main_edge_type = getMainEdgeType(edge_types[(u,v)])
        if main_edge_type is None:
            sys.stderr.write("WARNING: %s has no edge type. edge_types[%s]: %s. " % (str((u,v)),str((u,v)), str(edge_types[(u,v)])))
            sys.stderr.write("Evidence for edge: %s\n" % (str(evidence[(u,v)])))
            sys.stderr.write("\tSetting to 'physical'\n")
            #CHECK
            #sys.exit(1)
            # set this as the default for now
            main_edge_type = "physical"

        #if main_edge_type == '':
        #    raise NoEdgeTypeError("ERROR: %s,%s has no edge type. edge_types[%s,%s]: %s" % (u,v,u,v, str(edge_types[(u,v)])))

        if kwargs.get('no_k') is True:
            k_value = None
        else:
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
        edge_popup = buildEdgePopup(u,v, evidence, k=k_value)
        G.add_edge(gene_name_u,gene_name_v,directed=edge_dir[(u,v)],popup=edge_popup,k=k_value)

        attr_dict = {}
        if kwargs.get('case_study') is True:
            attr_dict['opacity'] = 0.8
        color = EDGE_COLORS[main_edge_type]
        edge_style = 'solid'
        edge_str = "%s-%s" % (u,v)
        if edge_str in graph_attr:
            if 'color' in graph_attr[edge_str]:
                color = graph_attr[edge_str]['color']
            if 'arrow_shape' in graph_attr[edge_str]:
                arrow_shape= graph_attr[edge_str]['arrow_shape']
            if 'edge_style' in graph_attr[edge_str]:
                edge_style = graph_attr[edge_str]['edge_style']

        G.add_edge_style(gene_name_u, gene_name_v, attr_dict=attr_dict,
                         directed=EDGE_DIR_MAP[main_edge_type], color=color, width=1.5, arrow_shape=arrow_shape,edge_style=edge_style)
    #G.set_data(metadata)
    #G.set_name(kwargs['graph_name'])
    return G


def get_color_luminance(color, hex_format=True):
    if hex_format:
        c = color.replace('#','')
        # this is from here: https://stackoverflow.com/a/29643643
        r,g,b = tuple(int(c[i:i+2], 16) for i in (0, 2, 4))
    else:
        print("ERROR: only hex_format is implemented for get_color_luminance()")
        sys.exit()
    # got this from here: https://stackoverflow.com/a/1855903
    luminance = (0.299 * r + 0.507 * g + 0.194 * b)/255.
    return luminance


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

    return main_edge_type


def buildNodePopup(n, pathswithnode=None, pathway_colors=None, chemical=None):
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
        # htmlstring += '<li><i>Recommended name</i>: %s</li>' % (db.get_description(uid))
        # get the alternate names from the AltName-Full namespace
        # alt_names = db.map_id(uid, 'uniprotkb', 'AltName-Full')
        # if len(alt_names) > 0:
        #     htmlstring += '<li><i>Alternate name(s)</i>: <ul><li>%s</li></ul></li>' %('</li><li>'.join(sorted(alt_names)))

        # #List EntrezGene IDs:
        # for e in db.map_id(uid, 'uniprotkb', 'GeneID'):
        #     htmlstring += '<br><b>Entrez ID:</b> <a style="color:blue" href="%s" target="EntrezGene">%s</a>' % (entrezurl%e, e)
#
#    #List IHOP link
#    ihopurl = 'http://www.ihop-net.org/UniPub/iHOP/?search=%s&field=UNIPROT__AC&ncbi_tax_id=9606' % (uid)
#    htmlstring += '<br><b>Search in </b> <a style="color:blue" href="%s" target="IHOP">IHOP</a>' % (ihopurl)

    if pathswithnode is not None:
        htmlstring += '<hr />'
        htmlstring += '<b>Paths</b>: %s<br>' %(','.join(str(i) for i in sorted(pathswithnode)))

    # if the pathways are specified for this node, add them to the list
    if pathway_colors is not None:
        htmlstring += "<hr /><b>Functions</b>:<ul>"
        for pathway in pathway_colors:
            if uid in pathway_colors[pathway]['prots']:
                pathway_link = '<a style="color:%s" href="%s">%s (%s)</a>' % (
                    pathway_colors[pathway]['color'], pathway_colors[pathway]['link'], pathway_colors[pathway]['name'], pathway)
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
def buildEdgePopup(t, h, evidence, k=None, family_ppi_evidence=None):
    #if t == "Q13315" and h == "Q01094":
    #    pdb.set_trace()
    annotation = '' 
    annotation+='<b>%s - %s</b></br>'%(','.join(sorted(uniprot_to_gene[t].split(','))), ','.join(sorted(uniprot_to_gene[h].split(','))))
    annotation+='<b>%s - %s</b></br>'%(t,h)
    if PPIWEIGHTS is not None:
        annotation+='<b>Weight</b>: %.3f</br>' % (PPIWEIGHTS[(t,h)])
    if k is not None:
        annotation+='<b>Edge Ranking</b>: %s' % (k)

    family_edge = True if len(t.split(',')) > 1 or len(h.split(',')) > 1 else False

    if family_ppi_evidence is not None and family_edge is True:
        annotation += '<hr /><h><b>Direct Sources of Evidence</b></h>'
        annotation += gs_utils.evidenceToHTML(t,h,evidence[(t,h)])
        if (t,h) in family_ppi_evidence:
            annotation += '<hr /><h><b>Sources of Evidence</b></h>'
            annotation += gs_utils.evidenceToHTML(t,h,family_ppi_evidence[(t,h)])
    else:
        annotation += '<hr /><h><b>Sources of Evidence</b></h>'
        annotation += gs_utils.evidenceToHTML(t,h,evidence[(t,h)])

    return annotation


# def post_to_graphspace(
#         version, interactome, rec_tfs_file, RESULTSPREFIX,
#         pathlinker_output_file, chemical, max_k=200, graphID='test',
#         postfix=None, tag=None, chemical_color_file=None):
#     # TODO change the call to post to graphspace so that the chemical map doesn't have to be read each time (slowing down posting)
#     # Could I update it so that the evidence doesn't have to be read each time as well?
#     # I would need to restructure post_to_graphspace_csbdb.py to be a class.
#     # new family network posting
#     command = "python src/graphspace/post_to_new_graphspace_evidence.py " + \
#               "--ppi %s " % (interactome) + \
#               "--sourcetarget %s " % (rec_tfs_file) + \
#               "--paths %s " % (pathlinker_output_file) + \
#               "--k-limit %d " % (max_k) + \
#               "--out-pref %s/graphspace/%s " % (RESULTSPREFIX, chemical) + \
#               "--username jeffl@vt.edu " + \
#               "--group ToxCast2 " + \
#               "--graph-name %s" % (graphID.replace(' ','_'))  # this should always be last because -casestudy will be appended to it if casestudy is true
#               #"--make-public " + \
#               # leave out the chemicalID for now since I need to update the post-to-graphspace script
#               #"--chemicalID %s " % (chemical) + \
#     # append the postfix to the graph name
#     if postfix:
#         command += postfix
#     # TODO add an option for case study graphs
#     casestudy = True
#     if casestudy:
#         # append '-casestudy' to the end of the graph name for now
#         # update: I'm using casestudy for the normal graph name, so leave it off for now
#         #command += '-casestudy'
#         command += " --casestudy "
#     else:
#         # apply the --human-rec-tfs option if this isn't a casestudy
#         command += " --human-rec-tfs $PWD/inputs/pathlinker-data/human-rec-tfs.txt "
#     ev_version = t_utils.get_ev_version(version)
#     command += " --version %s " % (ev_version)
#     # TEMP: The evidence file I used is different from the default:
#     command += "--evidence-file /data/jeff-law/svnrepo/data/interactions/compiled/2018_01/with-d2d/2018_01pathlinker-nokegg.tsv "
#     if tag:
#         command += " --tag %s " % (tag)
#     if chemical_color_file is not None:
#         command += " --graph-attr %s " % (chemical_color_file)
#         command += " --include-parent-nodes "
#     t_utils.runCommand(command, error_message="ERROR: Failed to post %s to graphspace\n"%rec_tfs_file)


def setup_parser():
    """
    """
    #parser = post_to_gs.setup_parser()
    ## Parse command line args.
    parser = argparse.ArgumentParser(description="Script to test for enrichment of the top predictions among given genesets")

    # general parameters
    group = parser.add_argument_group('Main Options')
    group.add_argument('--version',type=str,
            help="Version of the PPI to run. Options are: %s." % (', '.join(t_settings.ALLOWEDVERSIONS)))
    group.add_argument('--mapping-file', default="inputs/2018_01-toxcast-net/2018_01_uniprot_mapping.tsv",
            help='File to map to a different namespace. Network/edge IDs (uniprot ids) should be in the first column with the other namespace (gene name) in the second')
    group.add_argument('--revigo-file', 
            help="File containing the outputs of REVIGO for coloring the nodes in the graph")
    group.add_argument('--term-counts-file', 
            help="File containing the frequency of terms among the chemicals. Terms with a frequency < .75 will be selected.")
    group.add_argument('--ctd-support-file', 
            help="File containing the CTD phosphorylations interactions per chemical. Will be used to add a double border to nodes with support.")
    group.add_argument('--single-run','-S', type=str, action='append',
            help='Run only a single chemical. Can specify multiple chemicals by listing this option multiple times.')
    group.add_argument('--k-to-post', type=int, default=200,
            help='Value of k to test for significance. Multiple k values can be given.')
    group.add_argument('--forcepost',action='store_true', default=False,
            help='Force the network to be posted to graphspace even if json file already exists')   


    # posting options
    group = parser.add_argument_group('GraphSpace Options')
    group.add_argument('--username', '-U', type=str, 
                      help='GraphSpace account username to post graph to. Required')
    group.add_argument('--password', '-P', type=str,
                      help='Username\'s GraphSpace account password. Required')
    #group.add_argument('', '--graph-name', type=str, metavar='STR', default='test',
    #                  help='Graph name for posting to GraphSpace. Default = "test".')
    # group.add_argument('--out-pref', type=str, metavar='STR',
    #                   help='Prefix of name to place output files. ')
    group.add_argument('--name-postfix', type=str, default='',
                      help='Postfix of graph name to post to graphspace.')
    group.add_argument('--group', type=str,
                      help='Name of group to share the graph with.')
    group.add_argument('--case-study', action="store_true", default=False,
                      help='Use the Case study colors and labels (no k in labels, gray border around nodes)')
    group.add_argument('--make-public', action="store_true", default=False,
                      help='Option to make the uploaded graph public')
    # TODO implement and test this option
    #group.add_argument('--group-id', type=str, metavar='STR',
    #                  help='ID of the group. Could be useful to share a graph with a group that is not owned by the person posting')
    group.add_argument('--tag', dest="tags", type=str, action="append",
                      help='Tag to put on the graph. Can list multiple tags (for example --tag tag1 --tag tag2)')
    group.add_argument('--apply-layout', type=str,
                      help='Specify the name of a graph from which to apply a layout. Layout name specified by the --layout-name option. ' + 
                      'If left blank and the graph is being updated, it will attempt to apply the --layout-name layout.')
    group.add_argument('--layout-name', type=str, default='layout1',
                      help="Name of the layout of the graph specified by the --apply-layout option to apply. Default: 'layout1'")
    group.add_argument('--parent-nodes', action="store_true", default=False,
                      help='Use parent/group/compound nodes for the different node types')
    # group.add_argument('--graph-attr-file',
    #                    help='File used to specify graph attributes. Tab-delimited with columns: 1: style, 2: style attribute, ' + \
    #                    '3: nodes/edges to which styles will be applied separated by \'|\' (edges \'-\' separated), 4th: Description of style to add to node popup.')
    return parser


def parse_args(parser):
    # parse the command line arguments
    args = parser.parse_args()
    kwargs = vars(args)

    #for version in opts.version:
    if kwargs['version'] not in t_settings.ALLOWEDVERSIONS:
        print("ERROR: version '%s' not an allowed version. Options are: %s." % (kwargs['version'], ', '.join(t_settings.ALLOWEDVERSIONS)))
        sys.exit()

    return kwargs


if __name__ == '__main__':
    parser = setup_parser()
    kwargs = parse_args(parser)
    kwargs['chemicals'] = kwargs['single_run']
    main(**kwargs)


# def call_post_to_graphspace(chemicals, interactome, opts):
#     t_utils.checkDir("%s/graphspace" % (RESULTSPREFIX))
#     # write the color files of each chemical and return a dictionary of the chemical and its color file
#     # instead of writing all files here, write them one at a time so they can be posted in parallel with running MGSA
# #    chemical_color_files = {}
# #    if opts.mgsa_colors:
# #        chemical_color_files = mgsa_post_to_graphspace(chemicals, sig_cutoff=0.5, forced=opts.forcepost)
# #    elif opts.revigo_colors:
# #        chemical_color_files = revigo_post_to_graphspace(chemicals, forced=opts.forcepost)
#     # post everything to graphspace!
#     for chemical in tqdm(chemicals):
#         # get the chemical name. make sure it doesn't have any '%' or '[',']' as that will break posting
#         chemName = chemIDtoName[chemical].replace('%','').replace('[','').replace(']','')
#         rec_tfs_file = REC_TFS_FILE % (INPUTSPREFIX, chemical)
#         cyclinker_output_file = "%s/cyclinker/%s-paths.txt"%(RESULTSPREFIX, chemical)
#         output_json = "%s/graphspace/%s-graph.json" % (RESULTSPREFIX, chemical)
#         proteins, num_paths = t_utils.getProteins(paths=cyclinker_output_file, max_k=opts.k_to_post, ties=True)
#         if not opts.forcepost and os.path.isfile(output_json): 
#             print("%s already exists. Use --forcepost to overwrite it" % (output_json))
#         else:
#             # temporary fix to not re-post all chemical's networks
#             post_to_graphspace(interactome, rec_tfs_file, cyclinker_output_file, chemical, max_k=num_paths, graphID="%s_%s"%(chemName,chemical),
#                                postfix='-'+VERSION, tag=VERSION, chemical_color_file=opts.mgsa_colors)
#                                #postfix='-'+VERSION, tag=VERSION, chemical_color_file=chemical_color_files.get(chemical))


# def post_to_graphspace(interactome, rec_tfs_file, pathlinker_output_file, chemical, max_k=200, graphID='test', postfix=None, tag=None, chemical_color_file=False):
#     # TODO change the call to post to graphspace so that the chemical map doesn't have to be read each time (slowing down posting)
#     # Could I update it so that the evidence doesn't have to be read each time as well?
#     # I would need to restructure post_to_graphspace_csbdb.py to be a class.
#     # new family network posting
#     command = "python /home/jeffl/src/python/graphspace/trunk/graphspace-human/post_to_new_graphspace_evidence.py " + \
#               "--ppi %s " % (interactome) + \
#               "--sourcetarget %s " % (rec_tfs_file) + \
#               "--paths %s " % (pathlinker_output_file) + \
#               "--k-limit %d " % (max_k) + \
#               "--out-pref %s/graphspace/%s " % (RESULTSPREFIX, chemical) + \
#               "--username jeffl@vt.edu " + \
#               "--group ToxCast2 " + \
#               "--graph-name %s" % (graphID.replace(' ','_'))  # this should always be last because -casestudy will be appended to it if casestudy is true
#               #"--make-public " + \
#               # leave out the chemicalID for now since I need to update the post-to-graphspace script
#               #"--chemicalID %s " % (chemical) + \
#     # append the postfix to the graph name
#     if postfix:
#         command += postfix
#     # TODO add an option for case study graphs
#     casestudy = True
#     if casestudy:
#         # append '-casestudy' to the end of the graph name for now
#         # update: I'm using casestudy for the normal graph name, so leave it off for now
#         #command += '-casestudy'
#         command += " --casestudy "
#     else:
#         # apply the --human-rec-tfs option if this isn't a casestudy
#         command += " --human-rec-tfs $PWD/inputs/pathlinker-data/human-rec-tfs.txt "
#     ev_version = t_utils.get_ev_version(VERSION)
#     command += " --version %s " % (ev_version)
#     # TEMP: The evidence file I used is different from the default:
#     command += "--evidence-file /data/jeff-law/svnrepo/data/interactions/compiled/2018_01/with-d2d/2018_01pathlinker-nokegg.tsv "
#     if tag:
#         command += " --tag %s " % (tag)
#     if chemical_color_file is True:
#         chemical_color_files = mgsa_post_to_graphspace([chemical], sig_cutoff=0.5)
#         chemical_color_file = chemical_color_files[chemical] 
#         command += " --graph-attr %s " % (chemical_color_file)
#         command += "--include-parent-nodes "
#     t_utils.runCommand(command, error_message="ERROR: Failed to post %s to graphspace\n"%rec_tfs_file)


# def revigo_post_to_graphspace(chemicals, forced=False):
#     revigo_out_dir = "%s/stats/revigo" % (RESULTSPREFIX)

#     # define a set of colors to pull from
#     # TODO what if there are more than 12 pathways?
#     colors = t_settings.COLORS
#     out_dir = "%S/graphspace/colors" % (RESULTSPREFIX)
#     t_utils.checkDir(out_dir)
#     print("Writing REVIGO colors to %s for %d chemicals. (limit of %d colors)" % (out_dir, len(chemicals), len(colors)))
#     chemical_color_files = {}
#     for chemical in tqdm(sorted(chemicals)):
#         out_file = "%s/%s-colors.txt" % (out_dir, chemical)
#         chemical_color_files[chemical] = out_file
#         # if the file already exists and this isn't forced, then skip this round
#         if not forced and os.path.isfile(out_file):
#             continue

#         revigo_out_prefix = "%s/%s/%s" % (revigo_out_dir, chemical, chemical)
#         david_results_file = "%s-david-results.tsv" % (revigo_out_prefix) 
#         goterm_prots = {} 
#         goterm_pvals = {} 
#         goterm_to_id = {} 
#         # goterm + name are in the 2nd column, BF corrected p-value in the 11th
#         # annotated prots are in the 6th
#         for goterm_and_name, prots, pval in utils.readColumns(david_results_file, 2, 6, 11):
#             goterm_id, name = goterm_and_name.split("~")
#             goterm_prots[name] = prots.replace(', ', '|')
#             goterm_pvals[name] = pval
#             goterm_to_id[name] = goterm_id
#         revigo_selected_terms_file = "%s-selected-terms.txt" % (revigo_out_prefix)
#         revigo_selected_terms = utils.readItemList(revigo_selected_terms_file) 
#         #selected_terms = [goterm_to_id[name] for name in revigo_selected_terms]

#         goterm_popups = {}
#         link_template = "<a style=\"color:blue\" href=\"https://www.ebi.ac.uk/QuickGO/GTerm?id=%s\" target=\"DB\">%s</a>"
#         for goterm in revigo_selected_terms:
#             goterm_link = link_template % (goterm_to_id[goterm], goterm_to_id[goterm])
#             popup = "<b>QuickGO</b>:%s" % (goterm_link)
#             popup += "<br><b>p-value</b>: %0.2f" % (float(goterm_pvals[goterm]))
#             goterm_popups[goterm] = popup

#         write_colors_file(out_file, revigo_selected_terms, goterm_prots, goterm_popups)

#     return chemical_color_files


# def write_colors_file(out_file, functions, function_prots, function_popups=None, colors=None):
#     if colors is None:
#         colors = t_settings.COLORS
#     if len(functions) > len(colors):
#         print("\tWarning: # functions %d exceeds # colors %d. Limiting functions to the %d colors" % (len(functions), len(colors), len(colors)))
#         # limit the pathways to the number of colors
#         functions = functions[:len(colors)]

#     function_colors = {}
#     for i in range(len(functions)):
#         function_colors[functions[i]] = colors[i] 

#     # if no popup was provided, just show an empty popup
#     if function_popups is None:
#         function_popups = {function:"" for function in functions} 

#     print("\tWriting graphspace colors file: %s" % (out_file))
#     ## UPDATE 2017-07-13: create a compound or "parent" node for each goterm where the name of the parent node is the goterm name
#     ## TODO The 4th column is the description/popup of the parent node
#     # write the prots along with a color for posting to graphspace
#     with open(out_file, 'w') as out:
#         # first write the style of the individual nodes
#         out.write("#style\tstyle_val\tprots\tdescription\n")
#         out.write('\n'.join(['\t'.join(['color', function_colors[function], function_prots[function], '-']) for function in functions])+'\n')
#         # then write the parent nodes
#         out.write("#parent\tfunction\tprots\tpopup\n")
#         for function in functions:
#             out.write('\t'.join(['parent', function, function_prots[function], function_popups[function]])+'\n')
#         # then write the style of the parent node
#         out.write("#color\tcolor_val\tparent_node\tdescription\n")
#         out.write('\n'.join(['\t'.join(['color', function_colors[function], function, '-']) for function in functions])+'\n')


# def mgsa_post_to_graphspace(chemicals, sig_cutoff=0.5, forced=False):
#     mgsa_out_dir = "%s/stats/MGSA-out" % (RESULTSPREFIX)
#     #MsigDB  canonical pathways
#     #mgsa_input_file = "src/mgsa/c2.cp.v3.1.entrez.gmt"
#     mgsa_input_file = t_settings.MSIGDB_PATHWAYS
#     ##MsigDB  biological processes 
#     #mgsa_input_file = "src/mgsa/c5.bp.v5.0.entrez.gmt"

#     # define a set of colors to pull from
#     colors = t_settings.COLORS
#     print("Writing MGSA colors for %d chemicals. (limit of %d colors)" % (len(chemicals), len(colors)))
#     out_dir = "%s/graphspace/colors" % (RESULTSPREFIX)
#     t_utils.checkDir(out_dir)
#     # get the pathway_links from the mgsa input file
#     mgsa_pathway_nodes, pathway_links = mgsa_wrappers.read_mgsa_input_file(mgsa_input_file=mgsa_input_file)

#     chemical_color_files = {}
#     for chemical in tqdm(sorted(chemicals)):
#         out_file = "%s/%s-colors.txt" % (out_dir, chemical)
#         chemical_color_files[chemical] = out_file
#         # if the file already exists and this isn't forced, then skip this round
#         if not forced and os.path.isfile(out_file):
#             continue

#         mgsa_out_prefix = "%s/%s/%s" % (mgsa_out_dir, chemical, chemical)
#         mgsa_results_file = "%s-sig-results.txt" % (mgsa_out_prefix) 
#         # get the significant pathways sorted by the estimate 
#         sig_pathways = [(line[0], line[1][:4]) for line in utils.readColumns(mgsa_results_file, 1, 4)]

#         # get the nodes of that pathway in the cyclinker results 
#         pathway_prots_file = "%s-pathway-prots.txt" % (mgsa_out_prefix)
#         pathway_prots = {pathway: prots for pathway, prots in utils.readColumns(pathway_prots_file, 1, 2)}

#         pathway_popups = {}
#         for pathway, estimate in sig_pathways:
#             # TODO add these as a popup: estimate, pathway_links[pathway], 
#             pathway_link = "<a style=\"color:blue\" href=\"%s\" target=\"MsigDB\">%s</a>" % (pathway_links[pathway], pathway.replace("_", " "))
#             popup = "<b>Pathway</b>:%s" % (pathway_link)
#             popup += "<br><b>Posterior Probability</b>: %0.2f" % (float(estimate))
#             pathway_popups[pathway] = popup 

#         sig_pathways = [pathway for pathway, estimate in sig_pathways]

#         write_colors_file(out_file, sig_pathways, pathway_prots, pathway_popups)

#     return chemical_color_files

