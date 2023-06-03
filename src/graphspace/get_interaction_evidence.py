#!/usr/bin/python

# -----------------------------------------------
# -----------------------------------------------
import pdb
from optparse import OptionParser
import os
import sys
import networkx as nx
from src.utils import file_utils as utils
import glob
import re
from tqdm import tqdm
try:
    #sys.path.append('/home/jeffl/git-workspace/CSBDB')
    import query_csbdb
    import csbdb_utils
    db = csbdb_utils.db
except ImportError:
    sys.stderr.write("ImportError: Failed to import csbdb. " +
                     "This script imports csbdb_utils.py to map to different namespaces. " +
                     "You can download the CSBDB code using git:\n\tgit clone 'https://github.com/Murali-group/CSBDB.git'" +
                     "\nThen add CSBDB to your PYTHONPATH variable.\n")
    #raise
# parse_phosphosite.py is needed to parse the parsed BioPax phosphosite files. 
# currently it is located at src/python/DBParsers/trunk/PhosphositePlus/parse_phosphosite.py
try:
    import parse_phosphosite
except ImportError:
    sys.stderr.write("ImportError: This script retrieves phosphositeplus using the parse_phosphosite.py script located at src/python/DBParsers/trunk/PhosphositePlus/parse_phosphosite.py\n")
    #raise
# TODO only give an error if KEGG interactions are being used (in the evidence files? in the edges we're posting to graphspace?)
try:
    sys.path.append('/home/jeffl/svnrepo/src/python/DBParsers/trunk/KEGG')
    import kegg_to_graph
except ImportError:
    sys.stderr.write("ImportError: This script uses the kegg_to_graph.py script to get direction of KEGG edges " +
                     "(located at svnrepo/src/python/DBParsers/trunk/KEGG/kegg_to_graph.py).\nPlease add it to your PYTHONPATH variable\n")
    #raise
# -----------------------------------------------
# -----------------------------------------------

DATADIR = ''
VERSIONS = [
    '2018_01pathlinker',
    '2017_06pathlinker',
    '2017_01pathlinker',
    '2016_05pathlinker',
    '2015pathlinker',
]

# file which contains the mapping from the nodes (including family nodes) to the gene names
VERSIONMAPPINGFILE = {
    '2018_01pathlinker': "%s/interactions/compiled/2018_01/2018_01pathlinker-mapping.tsv",
    '2017_06pathlinker': "%s/interactions/compiled/2017_06/2017_06pathlinker-mapping.tsv",
    '2017_01pathlinker': "%s/interactions/compiled/2017_01/2017_01pathlinker-mapping.tsv",
    '2016_05pathlinker': "%s/interactions/compiled/2016_05/2016_05pathlinker-mapping.tsv",
    '2015pathlinker'   : "%s/interactions/compiled/2015/2015pathlinker-mapping.tsv",
}

VERSIONEVIDENCEFILE = {
    '2018_01pathlinker': "%s/interactions/compiled/2018_01/2018_01pathlinker.tsv",
    '2017_06pathlinker': "%s/interactions/compiled/2017_06/2017_06pathlinker.tsv",
    '2017_01pathlinker': "%s/interactions/compiled/2017_01/2017_01pathlinker.tsv",
    '2016_05pathlinker': "%s/interactions/compiled/2016_05/2016_05pathlinker.tsv",
    '2015pathlinker'   : "%s/interactions/compiled/2015/2015pathlinker.tsv",
}

# TODO allow the user to choose different sources from the same evidence file
# this would be more of a functionality for generate_human_interactome.py
VERSIONDATABASES = {
}

VERSIONSOURCES = {
    '2018_01pathlinker': {
        'CSBDB':'2018_01',  # the version of csbdb used to get the interactions
        'KEGG':'%s/interactions/kegg/2018-01-08/hsa/edge-files/', 
        'SPIKE':'%s/interactions/spike/edge-files-with-type/',       
        'PhosphositePlus': '/data/inputs/csbdb/data/2018_01/phosphosite/phosphositeplus',
    },   
    '2017_06pathlinker': {
        'CSBDB':'2017_06',  # the version of csbdb used to get the interactions
        'KEGG':'%s/interactions/kegg/2017-06-22/hsa/edge-files/', 
        'SPIKE':'%s/interactions/spike/edge-files/',       
        'PhosphositePlus': '/data/inputs/csbdb/data/2017_06/phosphosite/phosphositeplus',
    },   
    '2017_01pathlinker': {
        'CSBDB':'2017_01',  # the version of csbdb used to get the interactions
        'KEGG':'%s/interactions/kegg/2015-03-23/hsa/edge-files-with-family/', 
        'SPIKE':'%s/interactions/spike/edge-files/',       
        'PhosphositePlus': '/data/inputs/csbdb/data/2017_01/phosphosite/phosphositeplus',
    },   
    '2016_05pathlinker': {
        # use the file in SVN
        'CSBDB':'%s/csbdb/2016_05/2016-05-26-human_ppi-with-biogrid-r.txt',  
        # Does not have BioGrid in it (hence I left it as a build)
        #'CSBDB':'build_2016_05',  # the version of csbdb used to get the interactions. 
        #'NetPath':'%s/interactions/netpath/pathways/', NetPath is in CSBDB now
        'KEGG':'%s/interactions/kegg/2015-03-23/hsa/edge-files-with-family/', 
        'SPIKE':'%s/interactions/spike/edge-files/',       
        'PhosphositePlus' : '/data/inputs/csbdb/data/2016_05/phosphosite/parsed/phosphositeplus',
    },   
    '2015pathlinker': {
        # use the file in SVN
        'CSBDB':'%s/csbdb/2013-04-19-human_ppi-r.txt',
        'NetPath':'%s/interactions/netpath/pathways/',
        'KEGG':'%s/interactions/kegg/2015-03-23/hsa/edge-files/',
        'SPIKE':'%s/interactions/spike/edge-files/',       
    },   
}

NETPATHNAMES = {
    "NetPath_137" : "Advanced glycation end-products (AGE/RAGE)",
    "NetPath_1"   : "Alpha6 Beta4 Integrin",
    "NetPath_2"   : "Androgen receptor (AR)",
    "NetPath_12"  : "B cell receptor (BCR)",
    "NetPath_76"  : "Brain-derived neurotrophic factor (BDNF)",
    "NetPath_129" : "Corticotropin-releasing hormone (CRH)",
    "NetPath_4"   : "Epidermal growth factor receptor (EGFR)",
    "NetPath_25"  : "Follicle-stimulating hormone (FSH)",
    "NetPath_134" : "Fibroblast growth factor-1 (FGF1)",
    "NetPath_154" : "Gastrin",
    "NetPath_118" : "Ghrelin",
    "NetPath_10"  : "Hedgehog",
    "NetPath_5"   : "Inhibitor of differentiation (ID)",
    "NetPath_13"  : "Interleukin-1 (IL-1)",
    "NetPath_14"  : "Interleukin-2 (IL-2)",
    "NetPath_15"  : "Interleukin-3 (IL-3)",
    "NetPath_16"  : "Interleukin-4 (IL-4)",
    "NetPath_17"  : "Interleukin-5 (IL-5)",
    "NetPath_18"  : "Interleukin-6 (IL-6)",
    "NetPath_19"  : "Interleukin-7 (IL-7)",
    "NetPath_20"  : "Interleukin-9 (IL-9)",
    "NetPath_132" : "Interleukin-10 (IL-10)",
    "NetPath_147" : "Interleukin-11 (IL-11)",
    "NetPath_6"   : "Kit Receptor",
    "NetPath_22"  : "Leptin",
    "NetPath_3"   : "Notch",
    "NetPath_114" : "Oncostatin-M (OSM)",
    "NetPath_56"  : "Prolactin",
    "NetPath_21"  : "Receptor activator of nuclear factor kappa-B ligand (RANKL)",
    "NetPath_11"  : "T Cell Receptor (TCR)",
    "NetPath_7"   : "Transforming growth factor beta (TGF-beta) receptor",
    "NetPath_138" : "Tyrosine kinase with angiopoietins (TIE2/TEK)",
    "NetPath_9"   : "Tumor necrosis factor (TNF) alpha",
    "NetPath_23"  : "Thyroid-stimulating hormone (TSH)",
    "NetPath_24"  : "Thymic stromal lymphopoietin (TSLP)",
    "NetPath_26"  : "TNF-related weak inducer of apoptosis (TWEAK)",
    "NetPath_8"   : "Wnt",
}

def getEdgeDir(edges, evidence_file=None, split_family_nodes=False, add_ev_to_family_edges=False):
    edge_dir = {}
    print("Reading evidence file %s" % (evidence_file))
    evidence_lines = utils.readColumns(evidence_file, 1,2,3)

    # create a graph of the passed in edges 
    G = nx.Graph()
    G.add_edges_from(edges)
    if split_family_nodes or add_ev_to_family_edges:
        # a dictionary of the PPI evidence applied to the family edge is returned
        family_ppi_evidence = {}
        # also add the PPI edges in each family edge
        G.add_edges_from([(u,v) for U,V in edges for u in U.split(',') for v in V.split(',')])

        # get the family edges that are in the list of passed in edges. These can be ones that we created, so for now, keep them separate
        curr_family_edges = set([(U,V) for U,V in edges if len(U.split(',')) > 1 or len(V.split(',')) > 1])
        # get the family edges from the evidence file
        orig_family_edges = set()
        for line in evidence_lines:
            U = line[0]
            V = line[1]
            if len(U.split(',')) > 1 or len(V.split(',')) > 1:
                orig_family_edges.add((U,V))

        # only need to add the original family edges because only the original family edges will be in the evidence file
        G.add_edges_from(orig_family_edges)
        # initialize the dictionary
        for U,V in curr_family_edges.union(orig_family_edges):
            family_ppi_evidence[(U,V)] = {}
            family_ppi_evidence[(V,U)] = {}

        ppiCoveredByOrigFamily = mapPPItoFamily(G, orig_family_edges)
        ppiCoveredByCurrFamily = mapPPItoFamily(G, curr_family_edges)
        ppiCoveredByFamily = {}
        ppiCoveredByFamily.update(ppiCoveredByOrigFamily)
        ppiCoveredByFamily.update(ppiCoveredByCurrFamily)

    # initialize the dictionaries
    for t,h in G.edges():
        edge_dir[(t,h)] = False
        edge_dir[(h,t)] = False

    for t,h, directed in tqdm(evidence_lines):
        # header line of file
        #uniprot_a  uniprot_b   directed    interaction_type    detection_method    publication_id  source
        directed = True if directed == "True" else False

        # the set of edges for which to add this evidence to 
        edges_to_add_evidence = set([(t,h)])
        if split_family_nodes:
            for n1 in t.split(','):
                for n2 in h.split(','):
                    edges_to_add_evidence.add((n1,n2))
            # also check if only the head or tail was split
            #for n1 in t.split(','):
            #    edges_to_add_evidence.add((n1,h))
            #for n2 in h.split(','):
            #    edges_to_add_evidence.add((t,n2))

        for u,v in edges_to_add_evidence:
            # We only need to get the evidence for the edges passed in, so if this edge is not in the list of edges, skip it
            # G is a Graph so it handles undirected edges correctly
            if not G.has_edge(u,v):
                continue

            if directed is True:
                edge_dir[(u,v)] = True
            else:
                # if they are already directed, then set them as undirected
                if (u,v) not in edge_dir:
                    edge_dir[(u,v)] = False
                    edge_dir[(v,u)] = False

            if add_ev_to_family_edges and (u,v) in ppiCoveredByFamily:
                for U,V in ppiCoveredByFamily[(u,v)]:
                    # TODO this is an edge we created. The edges were created separately for directed and undirected. 
                    # Could be a problem if 4/5 edges were undirected and 1 was directed for this family edge
                    # for example: for a family node abcde, an undirected family edge abcde-x would be added from a-x, b-x, c-x, d-x
                    # if e->x was directed, there would not be a directed family edge abcde->x
                    # need to get that original decision
                    # for now, just set it as one of the PPI edges direction
                    if (U,V) not in edge_dir:
                        edge_dir[(U,V)] = directed

    return edge_dir


def getEvidence(edges, evidence_file=None, split_family_nodes=False, add_ev_to_family_edges=False):
    """
    *edges*: a set of edges for which to get the evidence for
    *split_family_nodes*: add evidence of family edges to pairwise interactions 
    *add_ev_to_family_edges*: add evidence of PPI edge covered by a family edge to the family edge. 
                              If this option is specified, then a dictionary of direct evidence along with
                              a dictionary of the ppi evidence applied to the family edge is returned
    returns a multi-level dictionary with the following structure
    edge: 
      db/source: 
        interaction_type: 
          detection_method: 
            publication / database ids
    for NetPath, KEGG and SPIKE, the detection method is the pathway name 
    NetPath, KEGG, and Phosphosite also have a pathway ID in the database/id set
    which follows convention of id_type:id (for example: kegg:hsa04611)
    """
    # dictionary of an edge mapped to the evidence for it
    evidence = {} 
    edge_types = {}
    edge_dir = {}

    if evidence_file is None:
        print("evidence_file not given. Returning empty sets")
        return evidence, edge_types, edge_dir
        # TODO use a default version or something

    print("Reading evidence file %s" % (evidence_file))
    evidence_lines = utils.readColumns(evidence_file, 1,2,3,4,5,6,7)

    # create a graph of the passed in edges 
    G = nx.Graph()
    G.add_edges_from(edges)
    if split_family_nodes or add_ev_to_family_edges:
        # a dictionary of the PPI evidence applied to the family edge is returned
        family_ppi_evidence = {}
        # also add the PPI edges in each family edge
        G.add_edges_from([(u,v) for U,V in edges for u in U.split(',') for v in V.split(',')])

        # get the family edges that are in the list of passed in edges. These can be ones that we created, so for now, keep them separate
        curr_family_edges = set([(U,V) for U,V in edges if len(U.split(',')) > 1 or len(V.split(',')) > 1])
        # get the family edges from the evidence file
        orig_family_edges = set()
        for line in evidence_lines:
            U = line[0]
            V = line[1]
            if len(U.split(',')) > 1 or len(V.split(',')) > 1:
                orig_family_edges.add((U,V))

        # only need to add the original family edges because only the original family edges will be in the evidence file
        G.add_edges_from(orig_family_edges)
        # initialize the dictionary
        for U,V in curr_family_edges.union(orig_family_edges):
            family_ppi_evidence[(U,V)] = {}
            family_ppi_evidence[(V,U)] = {}

        ppiCoveredByOrigFamily = mapPPItoFamily(G, orig_family_edges)
        ppiCoveredByCurrFamily = mapPPItoFamily(G, curr_family_edges)
        ppiCoveredByFamily = {}
        ppiCoveredByFamily.update(ppiCoveredByOrigFamily)
        ppiCoveredByFamily.update(ppiCoveredByCurrFamily)

    # initialize the dictionaries
    for t,h in G.edges():
        evidence[(t,h)] = {}
        evidence[(h,t)] = {}
        edge_types[(t,h)] = set()
        edge_types[(h,t)] = set()
        edge_dir[(t,h)] = False
        edge_dir[(h,t)] = False

    for t,h, directed, interactiontype, detectionmethod, pubid, source in tqdm(evidence_lines):
        # header line of file
        #uniprot_a  uniprot_b   directed    interaction_type    detection_method    publication_id  source
        directed = True if directed == "True" else False

        # the set of edges for which to add this evidence to 
        edges_to_add_evidence = set([(t,h)])
        if split_family_nodes:
            for n1 in t.split(','):
                for n2 in h.split(','):
                    edges_to_add_evidence.add((n1,n2))
            # also check if only the head or tail was split
            #for n1 in t.split(','):
            #    edges_to_add_evidence.add((n1,h))
            #for n2 in h.split(','):
            #    edges_to_add_evidence.add((t,n2))

        for u,v in edges_to_add_evidence:
            # We only need to get the evidence for the edges passed in, so if this edge is not in the list of edges, skip it
            # G is a Graph so it handles undirected edges correctly
            if not G.has_edge(u,v):
                continue

            evidence = addToEvidenceDict(evidence, (u,v), directed, source, interactiontype, detectionmethod, pubid)
            edge_types = addEdgeType(edge_types, (u,v), directed, source, interactiontype)
            if add_ev_to_family_edges and (u,v) in ppiCoveredByFamily:
                for U,V in ppiCoveredByFamily[(u,v)]:
                    # also add the evidence to the family edge
                    family_ppi_evidence = addToEvidenceDict(family_ppi_evidence, (U,V), directed, source, interactiontype, detectionmethod, pubid)
                    #if U == "P00533" and V == "P01111,P01112,P01116":
                    #    pdb.set_trace()
                    edge_types = addEdgeType(edge_types, (U,V), directed, source, interactiontype)
                    # TODO the edge was added separately for directed and undirected 
                    # need to get that original decision
                    # for now, just set it as one of the PPI edges direction
                    if (U,V) not in edge_dir:
                        edge_dir[(U,V)] = directed

            if directed is True:
                edge_dir[(u,v)] = True
            else:
                # if they are already directed, then set them as undirected
                if (u,v) not in edge_dir:
                    edge_dir[(u,v)] = False
                    edge_dir[(v,u)] = False

    if add_ev_to_family_edges:
        return evidence, family_ppi_evidence, edge_types, edge_dir
    else:
        return evidence, edge_types, edge_dir


def addEdgeType(edge_types, e, directed, source, interactiontype):
    if e not in edge_types:
        edge_types[e] = set()
    # add the edge type as well
    # direction was determined using the csbdb_interface.psimi_interaction_direction dictionary in CSBDB. 
    if not directed:
        # it would be awesome if we knew which edges are part of complex formation vs other physical interactions
        # is there an mi term for that?
        t,h = e
        edge_types[(t,h)].add('physical')
        if (h,t) not in edge_types:
            edge_types[(h,t)] = set()
        edge_types[(h,t)].add('physical')

    # include this to be able to reproduce evidence where there is no spike interaction type
    elif source == "SPIKE" and interactiontype == '':
        edge_types[e].add('spike_regulation')

    elif source == "KEGG" or source == "SPIKE":
        interactiontype = interactiontype.lower()
        if "phosphorylation" in interactiontype and "dephosphorylation" not in interactiontype:
            # Most of them are phosphorylation
            edge_types[e].add('phosphorylation')
        if "activation" in interactiontype:
            edge_types[e].add('activation')
        if "inhibition" in interactiontype:
            edge_types[e].add('inhibition')
        for edgetype in ["dephosphorylation", "ubiquitination", "methylation", "glycosylation"]:
            if edgetype in interactiontype:
                # TODO for now, just call the rest of the directed edges enzymatic. 
                edge_types[e].add('enzymatic')
    elif source == "Silverbush":
        edge_types[e].add('dir_predicted')
    else:
        # the rest is CSBDB
        #if "(phosphorylation reaction)" in interactiontype:
        # MI:0217 is safer because that will only ever match phosphorylation
        if "MI:0217" in interactiontype:
            # Most of the directed edges are phosphorylation
            edge_types[e].add('phosphorylation')
        else:
            # TODO for now, just call the rest of the directed edges enzymatic. 
            edge_types[e].add('enzymatic')

    return edge_types


def addToEvidenceDict(evidence, e, directed, source, interactiontype, detectionmethod, pubid):
    """ add the evidence of the edge to the evidence dictionary
    *pubids*: publication id to add to this edge. 
    """
    #if e not in evidence:
    #    evidence[e] = {}
    if source not in evidence[e]:
        evidence[e][source] = {}
    if interactiontype not in evidence[e][source]:
        evidence[e][source][interactiontype] = {}
    if detectionmethod not in evidence[e][source][interactiontype]:
        evidence[e][source][interactiontype][detectionmethod] = set()
    evidence[e][source][interactiontype][detectionmethod].add(pubid)

    if not directed:
        # add the evidence for both directions
        evidence = addToEvidenceDict(evidence, (e[1],e[0]), True, source, interactiontype, detectionmethod, pubid)

    return evidence


def mapPPItoFamily(G, family_edges):
    """
    returns a dictionary where the key is a PPI edge and the value is set the family edges that replaced it
    Note: a single PPI edge can be replaced by multiple family edges
    Consider the case xy->ab, xy->a, x->ab, and x->a. 
    x->a is replaced by xy->ab, xy->a, x->ab
    """
    # a mapping of a PPI edge to the family edge that is replacing it
    ppi_to_family = {}
    # first initialize the ppi_to_family dictionary
    for U,V in family_edges:
        for u in U.split(','):
            for v in V.split(','):
                # a single PPI edge can be replaced by multiple family edges
                ppi_to_family[(u,v)] = set()
                ppi_to_family[(v,u)] = set()
        for u in U.split(','):
            ppi_to_family[(u,V)] = set()
            ppi_to_family[(V,u)] = set()
        for v in V.split(','):
            ppi_to_family[(U,v)] = set()
            ppi_to_family[(v,U)] = set()

    # a mapping from the family edge to the set of PPI edges it is replacing
    #family_to_ppi = {}
    #for U,V in family_edges:
    #    family_to_ppi[(U,V)] = set()
    for U,V in family_edges:
        for u in U.split(','):
            for v in V.split(','):
                ppi_to_family[(u,v)].add((U,V)) 
                # TODO if dir doesn't trump undir, then this could be adding the directed evidence to an undirected edge
                # for example: consider the two edges v->u and u-v. if dir doesn't trump undir, then the v->u ev would be added to u-v
                # otherwise the u-v undirected edge wouldn't be in the interactome
                ppi_to_family[(v,u)].add((U,V)) 
        # add both direction so that undirected edges will also be able to map
        for u in U.split(','):
            ppi_to_family[(u,V)].add((U,V)) 
            ppi_to_family[(V,u)].add((U,V)) 
        for v in V.split(','):
            ppi_to_family[(U,v)].add((U,V)) 
            ppi_to_family[(v,U)].add((U,V)) 

    #return family_to_ppi, ppi_to_family
    return ppi_to_family


def parseArgs(args):
    # Parse command line args.
    usage = 'post_to_graphspace_human.py [options]\n'
    parser = OptionParser(usage=usage)
    parser.add_option('', '--datadir', type='str', metavar='STR', default='/home/jeffl/svnrepo/data',
                      help='Data directory from SVN. Required.')
    parser.add_option('','--version',type='string',default='2017_01pathlinker',
                      help='Version of interactome, which pulls from different directories.  Options are %s.  Default is "2017_01pathlinker".' % (','.join(VERSIONS)))
    # use the files specified by VERSIONEVIDENCEFILE
    parser.add_option('-o','--evidence-file',
                      help='File to write the sources of each human ppi to. If not specified, default is the default version: %s ' % VERSIONEVIDENCEFILE['2017_01pathlinker'])
    parser.add_option('','--split-family-nodes',action='store_true',default=False,
                      help='Split family nodes into protein nodes and add all pairwise combination PPIs of family edges.')
    parser.add_option('','--no-spike-edge-type',action='store_true',default=False,
                      help='Option included to reproduce old SPIKE files (before 2018-04-10) that did not have an interaction type column. ')
    parser.add_option('','--forced',action='store_true', default='False',
                      help="Force writing the evidence file")

    (opts, args) = parser.parse_args(args)

    if opts.version not in VERSIONS:
        sys.exit('ERROR: Version %s is not one of the possible versions. Possible versions are %s. Exiting.\n' % 
                 (opts.version,','.join(VERSIONS)))
    if opts.datadir is None:
        sys.exit('ERROR: data directory must be specified. Exiting.')

    return opts, args


def getEvidenceFile(version, datadir):
    print("Getting evidence file path")
    if version not in VERSIONS:
        sys.exit('ERROR: Version %s is not one of the possible versions. Possible versions are %s. Exiting.\n' % 
                 (version,','.join(VERSIONS)))
    return VERSIONEVIDENCEFILE[version] % datadir


def getMappingFile(version, datadir):
    print("Getting mapping file path")
    if version not in VERSIONS:
        sys.exit('ERROR: Version %s is not one of the possible versions. Possible versions are %s. Exiting.\n' % 
                 (version,','.join(VERSIONS)))
    return VERSIONMAPPINGFILE[version] % datadir


def main(version, datadir, evidence_file, forced=False, split_family_nodes=False, old_spike=False):
    global DATADIR
    DATADIR = datadir

    if not forced and os.path.isfile(evidence_file):
        print("Interaction evidence file %s already exists. Use --forced to overwrite it" % (evidence_file))
        return
    else:
        print("Compiling all interaction sources to the file: %s" % (evidence_file))
        evidence_file_dir = os.path.dirname(evidence_file) 
        if evidence_file_dir is not '' and not os.path.isdir(evidence_file_dir):
            print("Dir %s doesn't exist. Creating it" % (evidence_file_dir))
            os.makedirs(evidence_file_dir)

    sources = VERSIONSOURCES[version]
    # if we are getting the data from CSBDB, make sure schema is correct
    if 'CSBDB' in sources and '%' not in sources['CSBDB'] and \
            sources['CSBDB'] not in query_csbdb.SCHEMA_VERSIONS['all_versions']:
        sys.exit("ERROR: '%s' not an available schema version to download data from.\n\tAvailable versions: '%s'\n"%(sources['CSBDB'], "', '".join(query_csbdb.SCHEMA_VERSIONS['all_versions'])))
    # specify sources and pre-pend datadirectory.
    for s in sources:
        if '%' in sources[s]:
            sources[s] = sources[s] % (datadir)
    print("Using the following sources: %s" % ('; '.join(["%s: '%s'" % (s, sources[s]) for s in sources])))

    evidence, removed_evidence = compileEvidence(sources, split_family_nodes=split_family_nodes, old_spike=old_spike)

    removed_out_file = evidence_file + "removed-evidence.tsv"
    print("Writing evidence that was removed to: %s" % (removed_out_file))
    with open(removed_out_file, 'w') as out:
        # first write the header
        out.write("#%s\n" % ('\t'.join(['uniprot_a', 'uniprot_b', 'directed', 'interaction_type', 'detection_method', 'publication_id', 'source'])))
        # now write the evidence
        out_string = ""
        for (u,v) in tqdm(sorted(removed_evidence)):
            for ev in removed_evidence[(u,v)]:
                #out.write("\t".join([u, v, uniprot_to_gene[u], uniprot_to_gene[v]] + e))
                out_string += "\t".join([u, v] + [str(e) for e in ev]) + "\n"
        out.write(out_string)

    print("Writing evidence lines to %s" % (evidence_file))
    with open(evidence_file, 'w') as out:
        # first write the header
        out.write("#%s\n" % ('\t'.join(['uniprot_a', 'uniprot_b', 'directed', 'interaction_type', 'detection_method', 'publication_id', 'source'])))
        # now write the evidence
        out_string = ""
        for (u,v) in tqdm(sorted(evidence)):
            for ev in evidence[(u,v)]:
                #out.write("\t".join([u, v, uniprot_to_gene[u], uniprot_to_gene[v]] + e))
                out_string += "\t".join([u, v] + [str(e) for e in ev]) + "\n"
        out.write(out_string)

    # now map the uniprot IDs to gene names
    # TODO add an option to store the mapping file to csbdb_utils.py
    mapping_file = getMappingFile(version, datadir)
    if os.path.isfile(mapping_file):
        print("%s already exists. Not regenerating it" % (mapping_file))
        #uniprot_to_gene = utils.readDict(mapping_file, 1, 2)
    else:
        print("Mapping from uniprot to gene names")
        nodes = set([node for edge in evidence for node in edge])
        # split the family nodes into their proteins
        family_nodes = set()
        proteins = set()
        for n in nodes:
            if len(n.split(',')) > 1:
                family_nodes.add(n)
                proteins.update(set(n.split(',')))
            else:
                proteins.add(n)

        # map all of the proteins to gene names
        uniprot_to_gene = csbdb_utils.map_ids(proteins)

        # rejoin the gene names of the family nodes and add their mapping (i.e. uniprot_A,uniprot_B: gene_A,gene_B)
        print("%d family nodes" % (len(family_nodes)))
        for N in family_nodes:
            geneN = ','.join([uniprot_to_gene[n] for n in N.split(',')])
            uniprot_to_gene[N] = geneN

        # write the file so it can be used by other scripts (so the mapping doesn't have to be done by CSBDB each time which takes ~30 seconds for 17000 uniprot ids)
        print("Writing file mapping from uniprot to gene name: %s" % (mapping_file))
        with open(mapping_file, 'w') as out:
            out.write('\n'.join(["%s\t%s" % (n, uniprot_to_gene[n]) for n in uniprot_to_gene]))

    print("Finished")


def compileEvidence(evidence_files, split_family_nodes=False, old_spike=False):
    """
    gets the evidence for each interaction from each source in *evidence_files*
    returns a dictionary of each edge and the set of evidence tuples for it
    The tuples are as follows:
    directed, interaction_type, detection_method, publication_id, source
    """
    evidence = {}

    if 'SPIKE' in evidence_files:
        print("getting evidence for SPIKE")
        # TODO figure out a better way to distinguish between old way of parsing SPIKE and new way
        if old_spike is True:
            evidence = parseSPIKE_orig(evidence_files['SPIKE'], evidence)
        else:
            evidence = parseSPIKE(evidence_files['SPIKE'], evidence)

    if 'KEGG' in evidence_files:
        print("getting evidence for KEGG")
        evidence = parseKEGG(evidence_files['KEGG'], evidence, split_family_nodes=split_family_nodes)

    # make sure all of the edges are reviewed.
    # it's faster to access CSBDB as little as possible so checking the proteins after getting all of the edges is faster
    print("Ensuring all proteins in all edges are reviewed uniprot IDs")
    csbdb_utils.switch_to_current_schema()
    evidence, removed_evidence = removeNonReviewed(evidence)

    # TODO: fix the phosphositeplus ids in CSBDB 
    # no need to check the CSBDB proteins to see if they're reviewed because they came straight from the reviewed schema
    if 'CSBDB' in evidence_files:
        print("getting evidence for CSBDB")
        if 'PhosphositePlus' in evidence_files:
            evidence = parseCSBDB(evidence_files['CSBDB'], evidence, phosphositeplus_file=evidence_files['PhosphositePlus'])
        else:
            evidence = parseCSBDB(evidence_files['CSBDB'], evidence)

    return evidence, removed_evidence


def parse_psimi(mi_string):
    """ function to parse the PSI-MI string (such as 'psi-mi:"MI:0915"(physical association)')
    """
    mi = ''
    matchObj = re.search('(MI:\d+).*(\(.*\))',mi_string)
    if matchObj:
        mi = '%s %s' % (matchObj.group(1),matchObj.group(2))
    else:
        matchObj = re.search('(MI:\d+)',mi_string)
        if matchObj:
            mi = matchObj.group(1)

    return mi


def parseSPIKE(spike_dir, evidence):
    pathwayfiles = glob.glob('%s*-edges.txt' % (spike_dir))

    for f in tqdm(pathwayfiles):

        matchObj = re.search('%s(.*)_\d+_\d+.spike-edges.txt' % (spike_dir),f)
        if matchObj:
            pathwayname = matchObj.group(1)
            pathwayname = pathwayname.replace('_',' ')
        else:
            matchObj = re.search('%s(.*).spike-edges.txt' % (spike_dir),f)
            if matchObj:
                pathwayname = matchObj.group(1)
                pathwayname = pathwayname.replace('_',' ')
            else:
                matchObj = re.search('%s(.*)-edges.txt' % (spike_dir),f)
                if matchObj:
                    pathwayname = matchObj.group(1)
                    pathwayname = pathwayname.replace('_',' ')
                else:
                    sys.exit('ERROR! Cannot get pathway name from %s' % f)

        for t,h,interactiontype,integrity,pubmedids in utils.readColumns(f,1,2,5,6,7):
            # when the SPIKE files are parsed, regulations are directed, while interactions and complexes are undirected.
            if 'Regulation=' in interactiontype:
                directed = True
                e = (t,h)
            else:
                directed = False
                e = tuple(sorted((t,h)))
            # set the variables this way for now
            database = 'SPIKE'
            #interactiontype = ''
            detectionmethod = '%s (%s)' % (pathwayname, integrity)

            if e not in evidence:
                evidence[e] = set()

            for p in pubmedids.split(','):
                pubid = 'pubmed:%s' % p 
                evidence[e].add((str(directed), interactiontype, detectionmethod, pubid, database))

    return evidence


def parseSPIKE_orig(spike_dir, evidence):
    pathwayfiles = glob.glob('%s*-edges.txt' % (spike_dir))

    spikeG = nx.DiGraph()
    for f in pathwayfiles:
        lines = utils.readColumns(f,1,2)
        for t,h in lines:
            spikeG.add_edge(t,h)

    for f in tqdm(pathwayfiles):

        matchObj = re.search('%s(.*)_\d+_\d+.spike-edges.txt' % (spike_dir),f)
        if matchObj:
            pathwayname = matchObj.group(1)
            pathwayname = pathwayname.replace('_',' ')
        else:
            matchObj = re.search('%s(.*).spike-edges.txt' % (spike_dir),f)
            if matchObj:
                pathwayname = matchObj.group(1)
                pathwayname = pathwayname.replace('_',' ')
            else:
                matchObj = re.search('%s(.*)-edges.txt' % (spike_dir),f)
                if matchObj:
                    pathwayname = matchObj.group(1)
                    pathwayname = pathwayname.replace('_',' ')
                else:
                    sys.exit('ERROR! Cannot get pathway name from %s' % f)

        for t,h,integrity,pubmedids in utils.readColumns(f,1,2,5,6):
            # TODO when the SPIKE files are parsed, regulations are directed, while interactions and complexes are undirected.
            # regulations and interactions are indestinguishable in the edge files.
            # for now, we'll just have to use the fact if the reverse edge is in SPIKE or not.
            if spikeG.has_edge(h,t):
                directed = False
                e = tuple(sorted((t,h)))
            else:
                directed = True
                e = (t,h)
            # set the variables this way for now
            database = 'SPIKE'
            # TODO add an interactiontype. I think we could use physical association for undirected and direct interaction or something for directed
            interactiontype = ''
            detectionmethod = '%s (%s)' % (pathwayname, integrity)

            if e not in evidence:
                evidence[e] = set()

            for p in pubmedids.split(','):
                pubid = 'pubmed:%s' % p 
                evidence[e].add((str(directed), interactiontype, detectionmethod, pubid, database))

    return evidence


def parseKEGG(kegg_dir, evidence, split_family_nodes=False):
    global DATADIR
    print(DATADIR)
    keggnamefile = '%s/interactions/kegg/2015-03-23/hsa/HSA_PATHWAY_LIST_FORMATTED.txt' % (DATADIR)
    keggnames = {l[0]:l[1].replace('_',' ') for l in utils.readColumns(keggnamefile,1,2)}

    # links to KEGG pathway map
    kegg_map_link = 'http://www.kegg.jp/kegg-bin/show_pathway?'
    # links to KEGG pathway entry (evidence)
    kegg_entry_link = 'http://www.kegg.jp/dbget-bin/www_bget?pathway+'

    print("getting evidence for KEGG")
    if split_family_nodes:
        pathwaydir = '%s*-edges.txt' % (kegg_dir)
        pathwayfiles = glob.glob(pathwaydir)
    else:
        pathwaydir = '%s*-edges-family.txt' % (kegg_dir)
        pathwayfiles = glob.glob(pathwaydir)

    print("Getting edges from: %s" % (pathwaydir))
    #allEdges = set()
    for f in tqdm(pathwayfiles):
        if split_family_nodes:
            matchObj = re.search('%s(.*)-edges.txt' % (kegg_dir),f)
        else:
            matchObj = re.search('%s(.*)-edges-family.txt' % (kegg_dir),f)
        if matchObj:
            pathwayid = matchObj.group(1)
            pathwayname = keggnames.get(pathwayid,pathwayid)
        else:
            sys.exit('ERROR! Cannot get pathway name from %s' % f)
        for t,h,edgetype,edgesubtype in utils.readColumns(f,1,2,5,6):
            # use the kegg_to_graph.py function to get edge direction
            # it would be better if this information was stored in the edge files
            edgeDir = kegg_to_graph.determineEdgeDirection(t,h,edgetype,edgesubtype)
            # edgesubtype is a comma separated list of the following:
            # dir: activation, inhibition, phosphoryaltion, dephosphorylation, ubiquitination, methylation, glycosylation
            # indirect-effect, compound
            # undir: binding/association, dissociation, group-entry
            directed = False if edgeDir.lower() == 'undirected' else True

            e = (t,h)
            if not directed:
                e = tuple(sorted((t,h)))

            # edgetype seems to always be PPrel which isn't very helpful anyway. Just use the edgesubtype (for example: 'inhibition,phosphorylation')
            # also add that it is from KEGG
            interactiontype = edgesubtype
            #interactiontype = "KEGG %s" % (edgesubtype)
            detectionmethod = pathwayname
            keggid = "kegg:%s" % (pathwayid)
            database = 'KEGG'

            if e not in evidence:
                evidence[e] = set()

            evidence[e].add((str(directed), interactiontype, detectionmethod, keggid, database))

    return evidence


if __name__ == "__main__":
    opts, args = parseArgs(sys.argv)
    evidence_file = opts.evidence_file
    if evidence_file is None:
        evidence_file = getEvidenceFile(opts.version, opts.datadir)
    main(opts.version, opts.datadir, evidence_file, forced=opts.forced, 
            split_family_nodes=opts.split_family_nodes, old_spike=opts.no_spike_edge_type)
