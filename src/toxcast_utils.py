#!/usr/bin/python
""" File containing commonly used commands in ToxCast scripts
"""

import sys
import os
import re
import subprocess
import networkx as nx
import pandas as pd
from src import parse_toxcast_data as parse_toxcast_data
from src.utils import file_utils as utils


def copyToBaobabNodes(file_to_copy):
    import socket
    print("Copying %s to the 6 baobab nodes" % (file_to_copy))
    print("\tTip: If you are being asked for a password, setup ssh keys to allow passwordless ssh and scp (see the project ORG file)")
    ip_template = "192.168.200.%s"
    # use the localdisk storage on the baobab nodes to speed up read/write times
    copy_to_dir = re.sub("^/data", "/localdisk", os.path.dirname(os.path.abspath(file_to_copy))) 
    # loop through the 6 nodes
    for i in range(1,7):
        command = "ssh %s 'mkdir -p %s'; scp %s %s:%s" % (ip_template%i, copy_to_dir, os.path.abspath(file_to_copy), ip_template%i, copy_to_dir)
        if 'baobab' in socket.gethostname():
            subprocess.check_call(command, shell=True)
        else:
            runCommandOnBaobab(command)


def runCommandOnBaobab(command):
    # the -t command forces ssh to just try to run a command and then exit
    command = "ssh -t baobab.cbb.lan \"%s\"" % (command)
    print("Running: %s" % (command))
    subprocess.check_call(command, shell=True)


def setupEmailCommand(recipient, subject=None, message=None, baobab=False):
    """ Build a command for sending an email using an SMTP server with the 'sendEmail' package
    Currently very specific to my (Jeff's) setup
    """
    # default message if there were no errors
    if message is None:
        message = "Finished `date`"
    #print("Sending email to %s" % (recipient))
    # currently sent from wyatt
    # using my personal gmail (with a special app password)
    email_command = "sendemail -f <email> -t %s " % (recipient) + \
            "-u \"%s\" -m \"%s\" " % (subject, message) + \
            "-s smtp.gmail.com:587 -o tls=yes -xu <email> -xp \"<password>\""
    # if this command is going to be run on baobab, then ssh to wyatt first
    if baobab is True:
        email_command = "ssh -t jeffl@wyatt.cs.vt.edu \"%s\"" % (email_command)
    return email_command


def runCommand(command, show_output=True, quit=True, error_message=""):
    """ Small function to handle running bash subprocess
    *command*: The bash command to run on the system
    *show_output*: Default True. If False, output of the command will be hidden
    *quit*: Default True. If False, the script will not quit if the subprocess has an error
    *error_message*: Error message to show if the command failed
    """
    if show_output:
        print("Running:", command)
    try:
        if show_output:
            subprocess.check_call(command.split())
        else:
            subprocess.check_call(command.split(), stdout=open('/dev/null', 'w'), stderr=open('/dev/null', 'w'))
    except subprocess.CalledProcessError:
        if error_message:
            sys.stderr.write(error_message)
        if quit:
            raise


def loadToxcastData(interactome_file=None):
    """ load the toxcast data from parse_toxcast_data.py containing the hit matrix and various other mappings
    Also limit the data to the given interactome
    """
    print('loading chemical map')
    toxcast_data = parse_toxcast_data.ToxCastData()
    toxcast_data.load_data()
    #toxcast_data.get_chemical_types()
    if interactome_file is not None:
        # limit the data to the given interactome
        print("reading interactome from %s" % (interactome_file))
        G = build_network(interactome_file)
        print("\t%d nodes, %d edges" % (G.number_of_nodes(), G.number_of_edges()))
        toxcast_data.limit_to_interactome(G)
    #chemical_map.build_chemical_acc_map()
    return toxcast_data


def getToxcastHits(toxcast_data=None, limit_to_rec_tf=False):
    """ For each chemical, get a dictionary of the hit proteins (uniprot IDs).
    
    *toxcast_data*: ToxcastData object. If it doesn't exist, it will be created (i.e., files read)
    *limit_to_rec_tf*: limit the chemicals to those for which we have a signaling network
        (i.e., at least 1 hit rec and TF)
    """
    if toxcast_data is None:
        toxcast_data = loadToxcastData()
    hits_nonhits = toxcast_data.chemical_protein_hit
    if limit_to_rec_tf:
        hits_nonhits = {}
        # now limit the chemicals to those that have at least 1 rec and 1 TF hit
        for chem in toxcast_data.chemical_rec:
            hits_nonhits[chem] = toxcast_data.chemical_protein_hit[chem]
    # and get only the hit prots
    hits = {}
    for chem, prots in hits_nonhits.items():
        curr_hits = [] 
        for p, hit_val in prots.items():
            if hit_val == 1:
                curr_hits.append(p)
        hits[chem] = set(curr_hits)
    return hits


def get_ev_version(VERSION):
    """ Small hack to get the version of the evidence file used to build the data for this version of the interactome
    """
    if '2017_06' in VERSION:
        ev_version = "2017_06pathlinker"
    elif '2017_01' in VERSION:
        ev_version = "2017_01pathlinker"
    elif '2016_05' in VERSION:
        ev_version = "2016_05pathlinker"
    else:
        ev_version = "2018_01pathlinker"
    return ev_version


def build_network(net_file):
    ## Read the network from file
    net = nx.DiGraph()

    with open(net_file, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue
            items = [x.strip() for x in line.rstrip().split('\t')]
            id1 = items[0]
            id2 = items[1]
            eWeight = 1
            if(len(items) > 2):
                eWeight = float(items[2])
            net.add_edge(id1, id2, weight=eWeight)
    return net


def getChemicalTypeMaps(toxcast_data=None):
    """  Function to get the mapping from chemical DSS ID to chemical Type and vice versa
    """
    if toxcast_data is None:
        toxcast_data = parse_toxcast_data.ToxCastData()
    toxcast_data.get_chemical_types()
    return toxcast_data.chemDSStoTYPE, toxcast_data.chemTYPEtoDSS


def getChemicalNameMaps(toxcast_data=None):
    """  Function to get the mapping from chemical DSS ID to chemical Name and vice versa
    """
    if toxcast_data is None:
        toxcast_data = parse_toxcast_data.ToxCastData()
    toxcast_data.get_chemical_maps()
    #return chemical_map.chemDSStoName, chemical_map.chemNametoDSS
    return toxcast_data.chemIDtoName, toxcast_data.chemNametoID


def getChemicalCASMap(toxcast_data=None):
    """
    Map the chemical ID used in ToxCast to the CAS RN 
    """
    if toxcast_data is None:
        toxcast_data = parse_toxcast_data.ToxCastData()
    print("reading %s" % (toxcast_data.chemical_summary_file))
    df = pd.read_csv(toxcast_data.chemical_summary_file, header=0)
    cols = df.columns
    chemIDtoCAS  = dict(zip(df[cols[3]], df[cols[2]]))
    chemCAStoID  = dict(zip(df[cols[2]], df[cols[3]]))
    return chemIDtoCAS, chemCAStoID


def checkDir(directory):
    """ analagous to mkdir -p directory from the command line
    """
    if not os.path.isdir(directory):
        print("Dir %s doesn't exist. Creating it" % directory)
        try:
            os.makedirs(directory)
        except OSError:
            # if multiple parallel processes are running, they could be trying to create the same directory at the same time
            print("Dir %s already exists" % directory)


def getRecTFs(rec_tfs_file, costs=False):
    """ Get the receptors and TFs from a file. 
    uniprot_accession_number should be in the first column, with node_type (receptor or tf) in the second.
    """
    if costs is True:
        lines = utils.readColumns(rec_tfs_file,1,2,3,4)
        receptors = set([acc for acc, node_type, cost, zscore in lines if node_type == 'receptor'])
        tfs = set([acc for acc, node_type, cost, zscore in lines if node_type == 'tf'])
        costs = {acc:float(cost) for acc, node_type, cost, zscore in lines}
        zscores = {acc:float(zscore) for acc, node_type, cost, zscore in lines}
        return receptors, tfs, costs, zscores
    else:
        lines = utils.readColumns(rec_tfs_file,1,2)
        receptors = set([acc for acc, node_type in lines if node_type == 'receptor'])
        tfs = set([acc for acc, node_type in lines if node_type == 'tf'])
        return receptors, tfs


def build_chem_num_rectfs(chemicals, inputsprefix, rec_tfs_file_template="%s/rec-tfs/%s-rec-tfs.txt"):
    chem_num_rectfs = {}
    for chemical in chemicals:
        chem_rec, chem_tfs = getRecTFs(rec_tfs_file_template % (inputsprefix, chemical))
        chem_num_rectfs[chemical] = "%drec-%dtfs" % (len(chem_rec), len(chem_tfs))
    return chem_num_rectfs


def writeRecTFs(rec_tfs_file, receptors, tfs, costs=None, zscores=None):
    """ write the receptors and TFs to a file. 
    uniprot_accession_number should be in the first column, with node_type (receptor or tf) in the second.
    *weights*: dictionary of the cost for each rec/tf. Will be added as the third column if given.
               This cost can then be applied directly when running CycLinker as the cost of the super-source -> source and target -> super-target edges
    """
    with open(rec_tfs_file, 'w') as out:
        if costs is None:
            out.write("#node\tnode_type\n")
            # write the receptors
            out.write('\n'.join(["%s\treceptor" % (rec) for rec in receptors]) + '\n')
            # write the tfs
            out.write('\n'.join(["%s\ttf" % (tf) for tf in tfs]) + '\n')
        else:
            # also include the weight for each receptor and tf
            out.write("#node\tnode_type\tcost\tzscore\n")
            # write the receptors
            out.write('\n'.join(["%s\treceptor\t%0.6f\t%0.2f" % (rec, costs[rec], zscores[rec]) for rec in receptors]) + '\n')
            # write the tfs
            out.write('\n'.join(["%s\ttf\t%0.6f\t%0.2f" % (tf, costs[tf], zscores[tf]) for tf in tfs]) + '\n')


def getEdges(paths_file=None, paths=None, max_k=200, ties=False):
    """ Gets the edges from either the paths file, or the list of paths
    *returns*: a dictionary of the edge to the k value of the path it first came from
    """
    if paths is None:
        paths = getPaths(paths_file, max_k=max_k, ties=ties)
    edges = {}

    for k in range(len(paths)):
        path = paths[k]
        path = path.split('|')
        # get the edges from the paths file because the cyclinker edges file is not currently correct 
        for i in range(len(path)-1):
            edge = (path[i], path[i+1])
            if edge not in edges:
                # the paths are stored in order, so just use the index of the path +1 to get the actual k for the path
                edges[edge] = k+1

    return edges


def getPaths(paths_file, max_k=200, ties=False, scores=False):
    paths = []
    if scores:
        paths = {}
    last_score = None
    for k,score,path in utils.readColumns(paths_file, 1, 2, 3):
        # use 6 decimal places
        score = "%0.6f" % (float(score))
        if int(k) > max_k: 
            # if this path has the same path score as the previous path, then keep adding its proteins
            if ties and last_score == score:
                pass
            else:
                break 
        #path = path.split('|')
        if scores:
            paths[path] = score
        else:
            paths.append(path)
        last_score = score
    return paths


def getProteins(paths='', ranked_edges='', max_k=200, max_prots=None, ties=False):
    """ get the proteins of a network from the paths or ranked edges file
    The *ties* option  can only be used with the paths file
    The *ties* option will find the path score at the max_k, and then continue until that path score is passed
    *rec_tfs* option will include the receptors and TFs in the paths. Only works for the paths option
    """
    # get the proteins of the top k paths, or up to a certain number of proteins
    proteins = set()
    sources = set()
    targets = set()
    # keep track of the k or number of paths if ties is used
    num_paths = 0
    if paths:
        paths = getPaths(paths, max_k, ties)
        num_paths = len(paths) 
        for path in paths:
            path = path.split('|')
            # add the start and end of the path as source/target
            sources.add(path[0])
            targets.add(path[-1])
            # add each of the proteins in the path to the set of protiens
            proteins = proteins.union(set(path))
    else:
        # use the ranked edges file
        for p1,p2,k in utils.readColumns(ranked_edges, 1, 2, 3):
            if int(k) > max_k or len(proteins) > max_prots: 
                break 
            proteins.add(p1)
            proteins.add(p2) 

    # TODO add an option to also return family nodes
#    # remove the source/target family nodes from the set of proteins
#    for s in sources:
#        if len(s.split(',')) > 1:
#            proteins.remove(s)  
#    for t in targets:
#        if len(t.split(',')) > 1:
#            proteins.remove(t)  
#
#    # split family nodes into individual nodes
#    split_family_proteins = set()
#    for p in proteins:
#        split_family_proteins.update(set(p.split(',')))
#    proteins = split_family_proteins

    if ties:
        # return the total number of paths from keeping the ties
        return proteins, num_paths
    else:
        return proteins


def get_path_metrics(paths, max_k, ties=True, split_family_nodes=False):
    """
        loop through the paths to compute the number of paths each protein is in

    *split_family_nodes*: option to split family nodes in the counts 

    *returns* 
        1: a dictionary with uniprot IDs as the key, and counts as the value
        2: the set of edge tuples, 3: the set of receptors and tfs
    """
    num_paths_per_prot = {}
    # the number of times each protein shows up as an intermediate node
    num_inter_paths_per_prot = {}
    # the number paths each edge participates in
    num_paths_per_edge = {}
    last_score = 0
    rec_tfs = set()
    with open(paths, 'r') as paths:
        for line in paths:
            if line[0] == '#':
                continue
            line = line.rstrip().split('\t')
            k = int(line[0])
            score = "%0.6f" % (float(line[1]))
            path = line[2].split('|')
            if k > max_k:
                # if this path has the same path score as the previous path, then keep adding its proteins
                if ties and last_score == score:
                    pass
                else:
                    break
            last_score = score

            # keep track of the number of paths each protein is part of.
            for node in path:
                proteins = [node]
                # if this is a family node and the split_family_nodes option is set, then split it
                if split_family_nodes is True:
                    proteins = node.split(',')

                # if this node is not one of the sources or targets, then add it here.
                if path.index(node) > 0 and path.index(node) < len(path)-1:
                    for p in proteins:
                        if p not in num_inter_paths_per_prot:
                            num_inter_paths_per_prot[p] = 0 
                        num_inter_paths_per_prot[p] += 1
                for p in proteins:
                    if p not in num_paths_per_prot:
                        num_paths_per_prot[p] = 0 
                    num_paths_per_prot[p] += 1

            # add each edge to the set of edges
            for i in range(len(path)-1):
                edge = (path[i], path[i+1])
                if edge not in num_paths_per_edge:
                    num_paths_per_edge[edge] = 0 
                num_paths_per_edge[edge] += 1

            # add the first and last proteins in the path to the receptor and TF sets
            rec_tfs.add(path[0])
            rec_tfs.add(path[-1])

        return num_paths_per_prot, num_paths_per_edge, rec_tfs, num_inter_paths_per_prot


def computePval(chem_val, rand_vals, reverse=False):
    # TODO reverse should be the normal option. 
    if reverse:
        return sum([1 for i in range(0,len(rand_vals)) if rand_vals[i] >= chem_val]) / float(len(rand_vals)) 
    else:
        return sum([1 for i in range(0,len(rand_vals)) if rand_vals[i] <= chem_val]) / float(len(rand_vals)) 


def get_sig_chemicals(chemical_pvals, pval_col=5, sig_cutoff=0.05):
    """ *chemical_pvals* is the file output by compute_stat_sig.py
    """
    print("Getting the significant chemicals with a pval cutoff of %s from %s" % (str(sig_cutoff), chemical_pvals))
    sig_chemicals = []
    for chemical, xk_pval in utils.readColumns(chemical_pvals, 1, pval_col):
        if float(xk_pval) < sig_cutoff:
            sig_chemicals.append(chemical)
    print("%d chemicals are significant" % (len(sig_chemicals)))
    return sig_chemicals


def getPvals(resultsprefix, scope, sig_cutoff_type='FDR'):
    """ Function to retreive the pvalues for each chemical automatically.
        Currently only supports k 200 with FDR and BF pval corrections
    """
    print("Getting p-values for scope '%s'" % (scope))
    # get the significant chemicals
    # TODO add a k option and get the column automatically from the header line
    k = 200
    print("Using k %d" % (k))
    pvals_file = "%s/stats/stat-sig-%s/bfcorr_pval_qval.txt" % (resultsprefix, scope)
    with open(pvals_file, 'r') as file_handle:
        # example header line:
        #chemical   k10-BFcorr-pval k10-qval    k25-BFcorr-pval k25-qval    k50-BFcorr-pval k50-qval
        header = file_handle.readline().rstrip().split('\t')
    # TODO add an option to get the uncorrected pvals
    if sig_cutoff_type == 'BF':
        #pval_col = 10
        pval_col = header.index("k%d-BFcorr-pval" % k) + 1
    elif sig_cutoff_type == 'FDR':
        #pval_col = 11
        pval_col = header.index("k%d-qval" % k) + 1
    else:
        # TODO add the non-corrected p-value as an option
        print("please enter a valid value for --sig-cutoff-type. Valid options are: 'BF', 'FDR'")
        sys.exit(1)

    chemical_pvals = {chemical: float(k_pval) for chemical, k_pval in utils.readColumns(pvals_file, 1, pval_col)}
    return chemical_pvals


def getSigChemicals(resultsprefix, scope, sig_cutoff_type='FDR', sig_cutoff=0.05):
    """ Function to retreive the set of significant chemicals (response networks) automatically
    """
    chemical_pvals = getPvals(resultsprefix, scope, sig_cutoff_type='FDR')
    chemicals = [chemical for chemical in chemical_pvals if chemical_pvals[chemical] < sig_cutoff] 
    print("\t%d chemicals pass the significance cutoff of %0.2f with correction type of '%s'" % (len(chemicals), sig_cutoff, sig_cutoff_type))

    return chemicals


# Utils for working with family nodes
# -------------------------------------------------- 
def getFamilyNodes(version, interactome_file=None):
    """ Get the set of family nodes present in the interactome
    If the family-nodes.txt file already exists, then just read from it
    Otherwise read the interactome_file and write the family nodes to family-nodes.txt
    *returns*: the set of family nodes in the interactome
    """
    family_nodes_file = "inputs/%s/family-nodes.txt" % (version)
    if not os.path.isfile(family_nodes_file):
        print("%s does not exist. Getting the family nodes from the interactome" % (family_nodes_file))
        # get the set of family nodes from the interactome
        print("Reading the interactome from %s" % (interactome_file))
        family_nodes = set([N for U,V in utils.readColumns(interactome_file, 1, 2) for N in (U,V) if len(N.split(',')) > 1])
        print("Writing family nodes to %s" % (family_nodes_file))
        with open(family_nodes_file, 'w') as out:
            out.write('\n'.join(family_nodes))
    else:
        family_nodes = utils.readItemSet(family_nodes_file)

    return family_nodes


def getProteinToFamilyMapping(family_nodes):
    # mapping of a protein to the set of family nodes it's in
    proteinToFamily = {}
    for N in family_nodes:
        for n in N.split(','):
            if n not in proteinToFamily:
                proteinToFamily[n] = set()
            proteinToFamily[n].add(N)

    return proteinToFamily

