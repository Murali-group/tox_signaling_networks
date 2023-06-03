#! /usr/bin/env python

""" This script is for permuting a network and running edgelinker. This is meant to be run in parallel on a cluster (like baobab) and called from the toxcast master_script.py
"""

import sys, os
from optparse import OptionParser
#sys.path.append("/data/jeff-law/projects/2016-02-toxcast/src")
#sys.path.append("src")
#import toxcast_utils as t_utils
# make sure the base directory is the first thing in the path
if sys.path[0] != '':
    sys.path.insert(0, '')
from src.utils import file_utils as utils
from src import toxcast_utils as t_utils
from src import toxcast_settings as t_settings
from src import run_edgelinker
from src.stat_sig import compute_stat_sig
from tqdm import tqdm
# Python implementation of cycLinker
#import cycLinker
import networkx as nx
import permute_network
import random
import subprocess


def main():
    global VERSION, REC_TFS_PENALTY, EDGE_PENALTY
    opts = parse_args()

    if opts.version is not None:
        t_settings.set_version(opts.version)
        REC_TFS_PENALTY = t_settings.REC_TFS_PENALTY
        EDGE_PENALTY = t_settings.EDGE_PENALTY
        VERSION = opts.version
    else:
        REC_TFS_PENALTY = False
        EDGE_PENALTY = None
        VERSION = None

    # if a range or multiple indexes are specified, then call subprocesses of this script for each index.
    if opts.random_index_range or len(opts.random_index) > 1:
        # keep the --inputs-dir the relative path in the holmes directory for now
        command_template = "python -u src/permute/permute_and_edgelinker.py --interactome %s --inputs-dir %s --out-dir %s --k %d --random-index %%d --num-iterations %d " % \
                           (opts.interactome, opts.inputs_dir, opts.out_dir, opts.k, opts.num_iterations)
        # TODO make this automatic so I don't have to copy this each time I make an update
        if opts.write_score_counts:
            command_template += " --write-score-counts %s " % opts.write_score_counts
        if opts.undirected:
            command_template += " --undirected "
        if opts.swap_phys_sig_sep:
            command_template += " --swap-phys-sig-sep "
        if opts.split_by_weight:
            command_template += " --split-by-weight %s " % opts.split_by_weight
        if opts.forced:
            command_template += " --forced "
        if opts.cleanup:
            command_template += " --cleanup "
        if opts.version:
            command_template += " --version %s " % opts.version
        if opts.permute_rec_tfs:
            command_template += " --permute-rec-tfs %s " % opts.permute_rec_tfs
        #if opts.run_mgsa_random:
        #    command_template += " --run-mgsa-random "

        processes = set()
        if opts.random_index_range:
            indexes = range(opts.random_index_range[0],opts.random_index_range[1]+1)
        else:
            indexes = opts.random_index

        # run all of the specified indexes in parallel.
        for i in range(len(indexes)):
            index = indexes[i]
            # run these commands in parallel
            command = command_template % index
            print("Running", command)
            processes.add(subprocess.Popen(command.split()))
        # wait for them all to finish
        exit_codes = [p.wait() for p in processes]
        print("Finished." )
    else:
        # otherwise, run this index
        permute_and_run_edgelinker(opts, opts.random_index[0])


def permute_and_run_edgelinker(opts, random_index):
    if opts.write_score_counts:
        rand_scores_k = "%s/rand-networks/rand-%d-med-scores-k.txt" % (opts.write_score_counts, random_index)
        # if the final score counts file already exists, then don't do anything
        if os.path.isfile(rand_scores_k) and not opts.forced:
            print("%s already exists. Skipping." % (rand_scores_k))
            return
        chemical_k_scores = "%s/chemical-k-median-scores.txt" % (opts.write_score_counts)
        if not os.path.isfile(chemical_k_scores):
            print("Error: %s does not exist. Run compute_stat_sig.py with the --write-counts option to write it. Quitting" % (chemical_k_scores))
            return

    t_utils.checkDir("%s/networks" % (opts.out_dir))
    rec_tfs_file_template = "%s/rec-tfs/%%s-rec-tfs.txt" % (opts.inputs_dir) 
    chemicals = sorted(utils.readItemList("%s/chemicals.txt" % opts.inputs_dir, col=1))
    if opts.single_chem:
        chemicals = opts.single_chem

    if opts.permute_rec_tfs is not None:
        # if specified, "permute" the sets of receptors and tfs for each chemical instead of the interactome
        print("Writing random sets of rec/tfs for each chemical to %s" % (opts.out_dir))
        rec_tfs_file_template = "%s/%%s/%d-random-rec-tfs.txt" % (opts.out_dir, random_index)
        all_rec, all_tfs = t_utils.getRecTFs(opts.permute_rec_tfs)
        #chemical_num_rectfs_file = "%s/chemical_num_rectfs.txt" % (opts.inputs_dir)
        #lines = utils.readColumns(chemical_num_rectfs_file, 2, 3, 4)
        #for chem, num_rec, num_tfs in tqdm(lines):
        for chemical in tqdm(chemicals, disable=opts.verbose):
            out_file = rec_tfs_file_template % (chemical)
            if not os.path.isfile(out_file) or opts.forced:
                rec, tfs, costs, zscores = t_utils.getRecTFs(t_settings.REC_TFS_FILE % (opts.inputs_dir, chemical), costs=True)
                rec = list(rec)
                tfs = list(tfs)

                out_dir = "%s/%s" % (opts.out_dir, chemical)
                t_utils.checkDir(out_dir)
                random_rec = random.sample(all_rec, len(rec)) 
                # apply the costs to the random rec and tfs
                for i in range(len(rec)):
                    costs[random_rec[i]] = costs[rec[i]]
                    zscores[random_rec[i]] = zscores[rec[i]]
                random_tfs = random.sample(all_tfs, len(tfs)) 
                for i in range(len(tfs)):
                    costs[random_tfs[i]] = costs[tfs[i]]
                    zscores[random_tfs[i]] = zscores[tfs[i]]
                t_utils.writeRecTFs(out_file, random_rec, random_tfs, costs=costs, zscores=zscores)
        # use the original interactome
        permuted_network_out_file = opts.interactome
        print("Using the original interactome %s" % (permuted_network_out_file))
    else:
        # default is to permute the interactome
        permuted_network_out_file = '%s/networks/permuted-network%d.txt' % (opts.out_dir, random_index)
        if not os.path.isfile(permuted_network_out_file) or opts.forced:
            # don't log transform. The weights will be log transformed by the edgelinker code
            #G = cycLinker.readNetwork(opts.interactome, weight=True, logtransform=False)
            # UPDATE: 2017-12-07: try using the direction of the edges from the fourth column of the interactome instead of splitting based on if the edge is bidirected or not
            G = nx.DiGraph()
            dir_edges = []
            undir_edges = []
            lines = utils.readColumns(opts.interactome, 1, 2, 3, 4)
            if len(lines) == 0 :
                print("ERROR: interactome should have 4 columns: a, b, w, and True/False for directed/undirected. Quitting")
                sys.exit()
            for u, v, w, directed in lines:
                G.add_edge(u,v, weight=float(w))
                if directed.lower() in ["true", "t", "dir", 'directed']:
                    dir_edges.append((u,v))
                elif directed.lower() not in ["false", 'f', 'undir', 'undirected']:
                    print("ERROR: Unknown directed edge type '%s'. 4th column should be T/F to indicdate directed/undirected" % (directed.lower()))
                    print("Quitting.")
                    sys.exit()
                elif u < v:
                    undir_edges.append((u,v))

            if opts.undirected:
                # swap all edges as undirected edges
                permG = permute_network.permute_network(G.to_undirected(), num_iterations=opts.num_iterations)
                permG = permG.to_directed()
            elif opts.split_by_weight:
                # split the edges into bins by weight and swap the directed and undirected edges separately
                # if specified by the user
                permG = permute_network.permute_network(G, swap_phys_sig_sep=opts.swap_phys_sig_sep,
                                                        split_weight=opts.split_by_weight, num_iterations=opts.num_iterations)
            elif opts.swap_phys_sig_sep:
                # swap the directed and undirected edges separately
                permG = permute_network.permute_network(G, swap_phys_sig_sep=opts.swap_phys_sig_sep,
                                                        num_iterations=opts.num_iterations, edge_lists=(undir_edges, dir_edges))
            else:
                # if none of the options are specified, then swap everything as directed edges
                permG = permute_network.permute_network(G, num_iterations=opts.num_iterations)
            print("Writing %s" % (permuted_network_out_file))
            nx.write_weighted_edgelist(permG, permuted_network_out_file, comments='#', delimiter='\t')
        else:
            print("Using %s" % (permuted_network_out_file))

    # now run edgelinker on each of the chemicals using the permuted network
    # if version is netpath, use the different type of input file
    # TODO fix this 
    # PATHLINKERDATAVERSIONS
    #if 'kegg' in opts.inputs_dir or 'netpath' in opts.inputs_dir:
    #    rec_tfs_file_template = "%s/rec-tfs/%%s-nodes.txt" % (opts.inputs_dir) 
    in_files = []
    out_files = []
    for chemical in tqdm(chemicals, disable=opts.verbose):
        rec_tfs_file = rec_tfs_file_template % (chemical)
        in_files.append(os.path.abspath(rec_tfs_file))
        out_dir = "%s/%s" % (opts.out_dir, chemical)
        t_utils.checkDir(out_dir)
        out_pref = "%s/%d-random" % (out_dir, random_index)
        out_files.append(os.path.abspath(out_pref))
        # python implementation of edgelinker is taking too long. Switching to java for now.
        #run_write_edgelinker(permG, rec_tfs_file, opts.k, out_pref) 
        # run the java implementation of edgelinker below

    # write the in and out files to the networks dir
    edgelinker_in_files = '%s/networks/permuted-network%d-infiles.txt' % (opts.out_dir, random_index)
    with open(edgelinker_in_files, 'w') as out:
        out.write('\n'.join(in_files))
    edgelinker_out_files = '%s/networks/permuted-network%d-outfiles.txt' % (opts.out_dir, random_index)
    with open(edgelinker_out_files, 'w') as out:
        out.write('\n'.join(out_files))
    print("Running edgelinker on chemical %s: %s" % (chemical, out_pref))
    run_edgelinker.runEdgeLinker(permuted_network_out_file, edgelinker_in_files, edgelinker_out_files, opts.k, 
                 edge_penalty=EDGE_PENALTY, rec_tfs_penalty=REC_TFS_PENALTY, multi_run=True)

    if opts.write_score_counts:
        # now that edgelinker has been run on all of the chemical sources/targets,
        # get the path counts for the chemical network's path scores
        # import compute_stat_sig.py and run the code directly. This avoids the issues of re-importing the libraries from baobab
        print("Writing the counts for each of the scores for random index: '%d'" % (random_index))
        stat_sig = compute_stat_sig.StatSig(random_paths_dir=opts.out_dir, k_limit=opts.k, num_random=(random_index, random_index), out_dir=opts.write_score_counts)
        stat_sig.write_rand_counts(chemicals=chemicals, forced=opts.forced)
#        cmd = "python src/compute_stat_sig.py " + \
#              " --chemicals %s/chemicals.txt " % (opts.inputs_dir) + \
#              " --random-paths-dir %s/ " % (opts.out_dir) + \
#              " -P --k-limit %d " % (opts.k) + \
#              " --num-random %d %d" % (random_index, random_index) + \
#              " --group-by-prob " + \
#              " --write-rand-counts " + \
#              " --out-dir %s " % (opts.write_score_counts)
#        if opts.forced:
#            cmd += " --forced "
#        print(cmd)
#        subprocess.check_call(cmd.split())

    #if opts.run_mgsa_random:
    #    run_mgsa_random(random_index)

    if opts.cleanup:
        print("Deleting the generated permuted network and the edgelinker output files")
        if permuted_network_out_file != opts.interactome:
            os.remove(permuted_network_out_file)
        os.remove(edgelinker_in_files)
        # remove the individual output files
        for cyc_out_file in out_files:
            # # 2017-02-17 - temporarilly don't remove the paths file for running MGSA
            os.remove(cyc_out_file + "-paths.txt")
            os.remove(cyc_out_file + "-ranked-edges.txt")
        os.remove(edgelinker_out_files)


# # python implementation of edgelinker. 
# def run_write_edgelinker(G, annotations, k, out_pref):
#     # First make sure there is no 'source' or 'target' node from previous edgelinker runs
#     for node in ['source', 'target']:
#         if node in G.node:
#             G.remove_node(node)

#     ## Read in annotation file, identify sources and targets
#     sources, targets = cycLinker.readAnnotation(annotations)

#     ## Compute cycLinker Network
#     cycLinkerNetwork,paths = cycLinker.getCycLinkerNetwork(G, sources, targets, k, clip=False, logtransform=True)

#     ## write cycLinker results to file
#     cycLinker.printKSPGraph(out_pref+'-random-ranked-edges.txt',cycLinkerNetwork)
#     cycLinker.printKSPPaths(out_pref+'-random-paths.txt',paths)


def parse_args():
    parser = OptionParser()
    parser.add_option('','--version',type='string',
            help="Version of the PPI to run.  Options are: %s." % (', '.join(t_settings.ALLOWEDVERSIONS)))
    parser.add_option('','--interactome',type='string',default="/home/jeffl/svnrepo/data/interactomes/human/2016_05/2016-06-02-human-ppi-mod-sig-5_100-weighted.txt",
            help="Path to interactome to use. \n\tDefault = '/home/jeffl/svnrepo/data/interactomes/human/2016_05/2016-06-02-human-ppi-mod-sig-5_100-weighted.txt'")
    parser.add_option('','--inputs-dir', type='string',
                      help='Version dir containing the perturbed rec and tfs of each chemical. Used for --pval-hist option.')   
    parser.add_option('-o','--out-dir',type='string',default=".",
            help="Output directory to write the random networks to. \n\tDefault = '.'")   
    parser.add_option('','--write-score-counts',type='string',
            help="Output directory to write the path counts at specific scores. Will look for the 'chemical-k-median-scores.txt' in that dir and write the results to a 'rand-networks' dir.")   

    parser.add_option('-r','--random-index',type='int',action='append',
            help="Current permuted network index. If multiple indexes are specified, they will be run in parallel.")   
    parser.add_option('','--random-index-range',type='int',nargs=2,
            help="Range of indeces to permute network and run edgelinker. Will run a permute_and_cyclinker.py subprocess for each index in the specified range." +
                      "\n\tMeant to be run on a single node of cluster with the specified range equal to the # of cores on the node.")   

    parser.add_option('','--num-iterations',type='int',default='10',
            help="# of times to iterate through all edges and swap them. \n\tDefault = 10")   
    parser.add_option('','--k',type='int',default='200',
            help="# of paths to compute. \n\tDefault = 200")   
    parser.add_option('','--undirected',action='store_true',
            help="swap everything as undirected (default is to swap everything as directed edges)")   
    parser.add_option('','--swap-phys-sig-sep',action='store_true',
            help="swap the physical (undirected) and signaling (directed) edges seperately from each other")   
    parser.add_option('','--split-by-weight',type='float',
            help="swap the weights above and below the given weight separately")   
    parser.add_option('','--single-chem',action='append',
            help="specify one or more chemical to run on rather than the list of all of the chemicals")   
    parser.add_option('','--forced',action='store_true',
            help="force writing the interactome and running edgelinker even if the output files already exist.")   
    parser.add_option('','--cleanup',action='store_true',
            help="remove the generated permuted network after the edgelinker processes finish")   
    parser.add_option('','--run-mgsa-random',action='store_true',
            help="Run MGSA on the generated random network befor removing it.")   
    parser.add_option('','--permute-rec-tfs', action='store',
            help="Instead of permuting the interactome, select receptors and TFs uniformly at random from the file of receptors and TFs specified by this argument.")   
    parser.add_option('-v','--verbose', action='store_true', default=False,
            help="Option to display a progress bar when generating results. Default=False")   

    # parse the command line arguments
    (opts, args) = parser.parse_args()

    if opts.random_index_range and opts.random_index:
        print("Error: --random-index-range and --random-index cannot be used together")
        sys.exit(1)

    return opts


if __name__ == "__main__":
    main()
