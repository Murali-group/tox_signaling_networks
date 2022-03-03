#!/usr/bin/python
"""
File contains some functions related to breaking family nodes in an interactome, and add edge direction
"""

import sys
from tqdm import tqdm
import networkx as nx
import src.toxcast_utils as t_utils
import src.toxcast_settings as t_settings
import src.utils.file_utils as utils
# TODO
import get_interaction_evidence as getev


def getEvidenceVersion(version):
    # get the evidence from get_interaction_evidence.py
    if '2018_01' in version:
        evidence_version = "2018_01pathlinker"
    if '2017_06' in version:
        evidence_version = "2017_06pathlinker"
    elif '2017_01' in version:
        evidence_version = "2017_01pathlinker"
    else:
        evidence_version = "2016_05pathlinker"

    return evidence_version


def getUniprotToGeneMapping(version, datadir="/home/jeffl/svnrepo/data"):
    mapping_file = getev.getMappingFile(getEvidenceVersion(version), datadir)
    mapping = utils.readDict(mapping_file)

    return mapping


def addUndirPenalty(interactome_file, new_interactome_file, penalty=1.1):
    """ Increase the cost of undirected edges by decreasing the weight of those edges.
    *interactome_file*: tab-delimited with at least 4 columns: tail, head, weight, directed T/F
    *new_interactome_file*: path where new interactome will be written 

    The cost penalty would be adding log(penalty) to the cost of a path, 
    but it is simpler to apply that penalty directly to the edge weights by multiplying 1/penalty

    For example:
        w*(1/1.1) would decrease the weight by 10%, or increase the cost of a path by .041
    For an edge with weight 0.9:
        -log(0.9) = .046
        -log(0.9*1/1.1) = .087
    """
    print("Reading the interactome from %s" % (interactome_file))
    print("Adding a cost penalty to the undirected edges of log(%s)" % (penalty))
    print("\tApplying the penalty by multiplying each undirected edge by 1/%s (%0.3f)" % (
        penalty, 1/penalty))
    if penalty <= 1:
        sys.stderr.write("\nERROR, invalid penalty given. Must be >= 1. Quitting\n")
        sys.exit()
    print("Writing to %s" % (new_interactome_file))
    #t_utils.checkDir(os.path.dirname(new_interactome_file))
    out = ""
    with open(interactome_file, 'r') as f:
        # add the header line
        out += f.readline()
        for line in f:
            if line[0] == '#':
                continue
            line = line.rstrip().split('\t')
            if len(line) < 4 or line[3].lower() not in \
               ['true', 'false', 't','f', 'dir','undir']:
                sys.stderr.write("\nERROR: column indicating directed/undirected edge not found." +
                      "Must be 4th column. Quitting.\n")
                sys.exit()
            u,v,w,directed = line[:4]
            if directed.lower() in ['false', 'f', 'undir']:
                w = "%0.5e" % (float(w)*(1/penalty))
            out += "\t".join([u,v,w, directed])+'\n'
    with open(new_interactome_file, 'w') as f:
        f.write(out)
    return


## UPDATE 2017-06-26: If the sources/targets are contained within family nodes, then split those family nodes
def splitRecTFsFamilyNodes(chemicals, version, interactome_file):
    """
    """
    # leave some nodes as family nodes as that's how they are in the toxcast data
    map_family_to_prot = { 
        # FOS,JUN,FOSL1,FOSL2,JUNB,JUND,FOSB: FOS,JUN
        "P01100,P05412,P15407,P15408,P17275,P17535,P53539": ["P01100,P05412"],
        # FOS,JUN,SP1: FOS,JUN
        "P01100,P05412,P08047": ["P01100,P05412"],
        # FOS,JUN: FOS,JUN
        "P01100,P05412": ["P01100,P05412"],
        # TCF7,TCF7L1,TCF7L2,LEF1: TCF7,TCF7L1,TCF7L2,LEF1
        "P36402,Q9HCS4,Q9NQB0,Q9UJU2": ["P36402,Q9HCS4,Q9NQB0,Q9UJU2"],
        # FOXO3,FOXO4,FOXO1: FOXO3,FOXO4,FOXO1
        "O43524,P98177,Q12778": ["O43524,P98177,Q12778"],
    }

    rec_tfs_file = "inputs/%s/rec-tfs/%%s-rec-tfs.txt" % (version)
    interactomes_dir = "inputs/%s" % (version)
    t_utils.checkDir(interactomes_dir)
    new_interactome_file = "%s/%s-interactome.txt" % (interactomes_dir, version)
    # get the set of family nodes from the interactome
    print("Reading the interactome from %s" % (interactome_file))
    lines = utils.readColumns(interactome_file, 1,2,3)
    family_nodes = set([N for U,V,w in lines for N in (U,V) if len(N.split(',')) > 1])
    print("Splitting the source/target family nodes of all chemicals in the interactome and writing to %s" % (new_interactome_file))
    # set of family nodes to split from all chemicals
    family_to_split = {}
    for chemical in tqdm(chemicals):
        rec, tfs = t_utils.getRecTFs(rec_tfs_file % (chemical))
        for N in family_nodes:
            for n in rec.union(tfs):
                if n in N:
                    if N not in family_to_split:
                        family_to_split[N] = set()
                    family_to_split[N].add(n)
    # leave some tfs as family nodes because that's how they're listed in toxcast
    family_to_split.update(map_family_to_prot)

    split_rec = set()
    split_tfs = set()
    new_interactome = [] 
    all_new_edges = set()
    # it's a bit ad hoc because the weight of the family edge is the max of the individual edges, 
    # and now we're setting the edge weight of the split edges to be the max of the individual edges and the family edge
    new_edge_weights = {}
    #new_edge_ev = {} 
    # there could be multiple family edges contributing to a single edge
    for U,V,w in lines:
        new_edges = set()
        # split up the rec/tf family nodes 
        if U in family_to_split and V in family_to_split:
            split_rec.add(U) 
            split_tfs.add(V) 
            for u in family_to_split[U]:
                for v in family_to_split[V]:
                    new_edges.add((u,v))
        elif U in family_to_split:
            split_rec.add(U) 
            for u in family_to_split[U]:
                new_edges.add((u,V))
        elif V in family_to_split:
            split_tfs.add(V) 
            for v in family_to_split[V]:
                new_edges.add((U,v))
        # otherwise leave the edge as it is
        else:
            new_interactome.append((U,V,w))
            continue

        all_new_edges.update(new_edges)
        for (u,v) in new_edges:
            if (u,v) not in new_edge_weights:
                new_edge_weights[(u,v)] = set()
            new_edge_weights[(u,v)].add(float(w))
            # for now, don't write the evidence to each of the new networks to save on space
            # the evidence is present in the original interactome and the evidence file
            #if (u,v) not in new_edge_ev:
            #    new_edge_ev[(u,v)] = set()
            #new_edge_ev[(u,v)].update(set(ev.split('|')))

    for u,v in all_new_edges:
        w = max(new_edge_weights[(u,v)]) 
        #ev = '|'.join(new_edge_ev[(u,v)])
        new_interactome.append((u,v,"%0.6f"%w))

    # now write the new interactome
    print("Writing the new interactome with rec/tf family nodes split to %s" % (new_interactome_file))
    with open(new_interactome_file, 'w') as out:
        out.write('\n'.join(['\t'.join(line) for line in new_interactome])+'\n')

    # also write the family nodes that were split 
    mapping = getUniprotToGeneMapping(version)
    # also write the mapping from the rec/tf family node to the proteins it came from
    out_file = "inputs/%s/family-split-rec-tfs.txt" % (version)
    print("Writing a mapping of the split family rec/tfs and the protein hits they came from to: %s" % (out_file))
    with open(out_file, 'w') as out:
        out.write('\n'.join(["%s\t%s\t%s\t%s" % (N, '|'.join(family_to_split[N]), mapping[N], '|'.join([mapping[n] for n in family_to_split[N]])) for N in sorted(family_to_split)]) + '\n')

    print("A total of %d family nodes were split" % (len(family_to_split)))

    # add the zscore penalty to the few family nodes in the ToxCast data
    toxcast_family_nodes = [N[0] for N in map_family_to_prot.values()]
    addRecTFsFamilyNodes(chemicals, version, family_nodes=toxcast_family_nodes, costs=True)


def addEdgeDir(interactome_file, dir_trumps_undir=True, evidence_file=None, new_interactome_file=None):
    """ Add T/F for if the edge is directed or not as a 4th column 
    """
    if new_interactome_file is None:
        new_interactome_file = interactome_file
        print("Re-writing %s with edge direction as a fourth column" % (interactome_file))
    else:
        print("Reading %s and adding edge direction as a fourth column to %s" % (interactome_file, new_interactome_file))
    edges = set(utils.readColumns(interactome_file, 1, 2, 3))

    # ensure there are no self edges
    num_edges = len(edges)
    edges = [(u,v,w) for u,v,w in edges if u != v]
    if len(edges) != num_edges:
        print("%d self-edges removed" % (num_edges - len(edges)))

    edge_weights = {(u,v):w for u,v,w in edges}
    # also add the direction to the interactome
    print("Getting the edge direction of all edges in the interactome")
    if evidence_file is None:
        # try to get the edge direction automatically
        if '2018_01' in interactome_file:
            evidence_version = "2018_01pathlinker"
        if '2017_01' in interactome_file:
            evidence_version = "2017_01pathlinker"
        else:
            evidence_version = "2016_05pathlinker"
        evidence_file = getev.getEvidenceFile(evidence_version, t_settings.DATADIR)
    # the evidence file has if the edges are directed or not. 
    # get that information here using get_interaction_evidence.py
    edge_dir = getev.getEdgeDir(edge_weights.keys(), evidence_file, split_family_nodes=True, add_ev_to_family_edges=True)

    # UPDATE 2017-12-07: After splitting the family receptors and TFs,
    # some of the undirected edges that are trumped by directed edges are back again. Make sure they're removed
    G = nx.Graph()
    dirG = nx.DiGraph()
    for u,v in edge_weights:
        if edge_dir[(u,v)] is True:
            dirG.add_edge(u,v)
        else:
            G.add_edge(u,v)

    # now remove trumped edges
    trumped_edges = 0
    for u,v in dirG.edges():
        if G.has_edge(u,v):
            trumped_edges += 1
            G.remove_edge(u,v)

    print("%d undir edges trumped by dir edges" % (trumped_edges))
    print("%d directed, %d undirected edges" % (dirG.number_of_edges(), G.number_of_edges()))
    combG = nx.DiGraph()
    combG.add_edges_from(dirG.edges(), type="Dir")
    combG.add_edges_from(G.to_directed().edges(), type="Undir")
    print("%d total edges" % (combG.number_of_edges()))

    #new_edges = dirG.edges()
    #for u,v in G.edges():
    #    new_edges.append((u,v))
    #    new_edges.append((v,u))
    #if len(new_edges) != len(set(new_edges)):
    #    print("ERROR: there are duplicates")
    if combG.number_of_edges() != (dirG.number_of_edges() + (G.number_of_edges()*2)):
        print("ERROR: len(new_edges) != dirG.number_of_edges() + (G.number_of_edges()*2):")
        print("%d != %d + (%d*2)" % (combG.number_of_edges(), dirG.number_of_edges(), G.number_of_edges()))
        sys.exit()

    # now write the interactome file again
    with open(new_interactome_file, 'w') as out:
        for u,v in combG.edges():
            #dir_string = "Dir" if edge_dir[(u,v)] else "Undir"
            w = edge_weights[(u,v)] if (u,v) in edge_weights else edge_weights[(v,u)]
            out.write("%s\t%s\t%s\t%s\n" % (u,v,w,combG.edge[u][v]['type']))


def addRecTFsFamilyNodes(chemicals, version, interactome_file=None, family_nodes=None, costs=True):
    """ add family nodes to the set of sources and targets for each chemical
    for now, just add every family node that has a hit protein in it
    *interactome_file* used to get the set of family nodes in the interactome. not used if family_nodes are passed in
    *family_nodes* family nodes to check for contained rec/tfs hits. If None, all family nodes in the interactome will be checked
    for example: inputs/version/rec-tfs/%s-rec-tfs.txt
    """
    rec_tfs_file = "inputs/%s/rec-tfs/%%s-rec-tfs.txt" % (version)

    if family_nodes is None and interactome_file is not None:
        family_nodes = t_utils.getFamilyNodes(version, interactome_file)
    # mapping of a protein to the set of family nodes it's in
    proteinToFamily = t_utils.getProteinToFamilyMapping(family_nodes)
    familyToProtein = {}
    print("Updating the set of sources and targets with family nodes for %d chemicals" % (len(chemicals)))
    rec_added = set()
    tfs_added = set()
    rec_tf_costs = {}
    rec_tf_zscores = {}
    for chemical in tqdm(chemicals):
        if costs is True:
            rec, tfs, rec_tf_costs, rec_tf_zscores = t_utils.getRecTFs(rec_tfs_file % (chemical), costs=True)
        else:
            rec, tfs = t_utils.getRecTFs(rec_tfs_file % (chemical))
            # just build empty dictionaries for the cost and zscore dictionaries and don't write them to the file
            rec_tf_costs = {acc: '' for acc in rec.union(tfs)}
            rec_tf_zscores = {acc: '' for acc in rec.union(tfs)}

        new_rec = set()
        new_tfs = set()
        for acc in rec.union(tfs):
            if acc in proteinToFamily:
                # add the family nodes that have this receptor/tf protein in it
                if acc in rec:
                    new_rec.update(proteinToFamily[acc])
                    for N in proteinToFamily[acc]:
                        if N not in familyToProtein:
                            familyToProtein[N] = set() 
                        familyToProtein[N].add(acc)
                else:
                    new_tfs.update(proteinToFamily[acc])
                # Add the cost of the protein to the family node
                for N in proteinToFamily[acc]:
                    rec_tf_costs[N] = float(rec_tf_costs[acc])
                    rec_tf_zscores[N] = float(rec_tf_zscores[acc])
                    for N in proteinToFamily[acc]:
                        if N not in familyToProtein:
                            familyToProtein[N] = set() 
                        familyToProtein[N].add(acc)
        # keep track of the ones added 
        rec_added.update(new_rec)
        tfs_added.update(new_tfs)

        if len(rec.union(new_rec)) == len(rec) and len(tfs.union(new_tfs)) == len(tfs):
            # skip writing the file if it's the same
            continue
        if costs is True:
            t_utils.writeRecTFs(rec_tfs_file % (chemical), rec.union(new_rec), tfs.union(new_tfs), costs=rec_tf_costs, zscores=rec_tf_zscores)
        else:
            t_utils.writeRecTFs(rec_tfs_file % (chemical), rec.union(new_rec), tfs.union(new_tfs))

    mapping = getUniprotToGeneMapping(version)
    # also write the mapping from the rec/tf family node to the proteins it came from
    out_file = "inputs/%s/family-mapping-rec-tfs.txt" % (version)
    print("Writing a mapping of the new family rec/tfs and the protein hits they came from to: %s" % (out_file))
    with open(out_file, 'w') as out:
        out.write('\n'.join(["%s\t%s\t%s\t%s" % (N, '|'.join(familyToProtein[N]), mapping[N], '|'.join([mapping[n] for n in familyToProtein[N]])) for N in sorted(familyToProtein)]))

    print("A total of %d family rec and %d family tfs were added to the chemicals" % (len(rec_added), len(tfs_added)))
# -------------------------------------------------- 
