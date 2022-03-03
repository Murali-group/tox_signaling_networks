#!/usr/bin/env python

## Jeff Law, March 2022
## Master Script for running Toxcast

print("Importing libraries")
from optparse import OptionParser,OptionGroup
from src.utils import file_utils as utils
from src import toxcast_utils as t_utils
from src import interactome_utils as i_utils
# script to hold the different versions and settings and such
from src import toxcast_settings as t_settings
# used for running edgelinker in a parallel manner
from src import run_edgelinker
import os
import sys
import subprocess
import re
#import src.map_chemicals as map_chemicals
# Add /data/jeff-law/tools/tqdm to your python path to give the progress bar functionality.
from tqdm import tqdm
# for checking to see if this script is being run on baobab
import socket


TOXCAST_DATA = None
chemIDtoName, chemNametoID = t_utils.getChemicalNameMaps()


def main(args):
    global VERSION, DATADIR, INTERACTOMES
    global RESULTSPREFIX, INPUTSPREFIX, REC_TFS_FILE
    global PATHLINKERDATAVERSIONS, TOXCAST_DATA
    global REC_TFS_PENALTY
    ## parse arguments
    opts = parseArguments(args)

    for version in opts.version:
        VERSION = version

        # settings are now stored in the src/toxcast_settings.py script
        INPUTSPREFIX, RESULTSPREFIX, interactome = t_settings.set_version(VERSION)
        REC_TFS_PENALTY = t_settings.REC_TFS_PENALTY
        EDGE_PENALTY = t_settings.EDGE_PENALTY
        DATADIR = t_settings.DATADIR
        #PATHLINKERDATAVERSIONS = t_settings.PATHLINKERDATAVERSIONS 
        REC_TFS_FILE = t_settings.REC_TFS_FILE 

        # create the directories if they don't exist
        t_utils.checkDir(RESULTSPREFIX)
        t_utils.checkDir(INPUTSPREFIX)

        # check if the inputs are there
        if opts.forceinputs or not os.path.isdir("%s/rec-tfs" % INPUTSPREFIX):
            # load the toxcast data
            TOXCAST_DATA = t_utils.loadToxcastData(t_settings.INTERACTOMES[VERSION])
            if VERSION in t_settings.ZSCORE_PENALTY:
                print("Adding a zscore cost to the sources/targets")
                include_zscore_weight = True 
            else:
                include_zscore_weight = False 
            # generate the perturbed rec and tf files as well as the chemicals.txt file
            TOXCAST_DATA.write_chemical_perturbed_rec_tfs(
                "%s/chemicals.txt" % INPUTSPREFIX, "%s/rec-tfs" % INPUTSPREFIX,
                include_zscore_weight=include_zscore_weight)
        else:
            TOXCAST_DATA = t_utils.loadToxcastData(t_settings.INTERACTOMES[VERSION])
        # these chemicals have at least 1 hit rec and tf
        chemicals = TOXCAST_DATA.chemical_rec.keys()

        # If we're penalizing the undirected edges, add that penalty now
        if VERSION in t_settings.UNDIR_PENALTIES and not os.path.isfile(interactome):
            i_utils.addUndirPenalty(
                t_settings.INTERACTOMES[VERSION], interactome, penalty=t_settings.UNDIR_PENALTIES[VERSION])

        # kind of like a test run before starting everything
        if opts.single_run:
            #chemicals = [chemicals[0]]
            #chemicals = [chemNametoID['Tiratricol']]
            chemicals = sorted(opts.single_run)
            # if the user specified to use the 'single_chem.txt' file, use that
            if 'file' in opts.single_run:
                del chemicals[chemicals.index('file')]
                single_chem_file = "%s/single_chem.txt" % INPUTSPREFIX
                if os.path.isfile(single_chem_file):
                    print("Using %d chemicals plus those listed in '%s'" % (len(chemicals), single_chem_file))
                    chemicals += sorted(utils.readItemList(single_chem_file))
                else:
                    print("file '%s' not found." % (single_chem_file))

        print("Using version '%s' with %d chemicals" % (VERSION, len(chemicals)))

        if opts.response_networks:
            generateResponseNetworks(chemicals, interactome, max_k=opts.edgelinker_k, pathlinker=opts.pathlinker, 
                                     edge_penalty=EDGE_PENALTY, rec_tfs_penalty=REC_TFS_PENALTY, forced=opts.forcealg)

        # if opts.graphspace:
        #     call_post_to_graphspace(chemicals, interactome, opts) 

        scopes = opts.scope
        if scopes is None:
            scopes = [scopes]
        for scope in scopes:
            if opts.random_networks:
                if 'permute' in scope:
                    print("Permuting the interactome %d times, and then running edgelinker on each of the permuted networks" % (10))
                else:
                    print("Permuting the set of receptors and TFs for each chemical, then running edgelinker on the original interactome")
                    print("Using scope '%s'" % (scope))
                edgelinkerPermuteNetworks(interactome, scope, opts.num_random_networks, opts.edgelinker_k, opts.k_to_test, opts.forcealg,
                                        opts.parallel_processes, opts.super_computer, num_iterations=10,
                                        printonly=opts.printonly, cleanup=opts.cleanup, single_chem=opts.single_run,
                                        send_email=opts.send_email)

            # TODO update this option for the new statistical test
            if opts.stat_sig:
                print("Using scope '%s'" % (scope))
                compute_stat_sig(scope, opts.num_random_networks, opts.k_to_test, forced=opts.forcesig)

            if opts.gen_figs:
                gen_figs(scope, opts.k_to_test, forced=opts.forcesig)

            # keep track of this for later
            total_num_chemicals = len(chemicals)
            curr_chemicals = chemicals

            # get only the chemicals that have a significant response network
            if opts.sig_cutoff:
                sig_chemicals = t_utils.getSigChemicals(RESULTSPREFIX, scope, sig_cutoff_type=opts.sig_cutoff_type, sig_cutoff=opts.sig_cutoff)
                nonsig_chemicals = set(chemicals).difference(sig_chemicals)
                # write these to a file
                sig_chemicals_file = "%s/sig-chemicals.txt" % (INPUTSPREFIX)
                print("Writing the current set of significant chemicals to %s" % (sig_chemicals_file))
                with open(sig_chemicals_file, 'w') as out:
                    header = "# chemicals from scope: '%s', sig_cutoff_type: '%s', sig_cutoff: '%s'\n" % (scope, opts.sig_cutoff_type, opts.sig_cutoff)
                    out.write(header)
                    out.write('\n'.join(sig_chemicals))
                nonsig_chemicals_file = "%s/nonsig-chemicals.txt" % (INPUTSPREFIX)
                print("Writing the current set of non-significant chemicals to %s" % (nonsig_chemicals_file))
                with open(nonsig_chemicals_file, 'w') as out:
                    header = "# chemicals from scope: '%s', sig_cutoff_type: '%s', sig_cutoff: '%s'\n" % (scope, opts.sig_cutoff_type, opts.sig_cutoff)
                    out.write(header)
                    out.write('\n'.join(nonsig_chemicals))
                # use the sig chemicals as the set of chemicals for the rest of the scripts
                curr_chemicals = sig_chemicals

            # add a prefix to the output files showing which significance settings were used for this set of chemicals
            sig_pref = 'all'
            if opts.sig_cutoff:
                sig_pref = "%s-%s-%s" % (scope, opts.sig_cutoff_type, str(opts.sig_cutoff).replace('.','_'))

            if opts.overlap:
                # TODO allow generating both with a single call to the master script
                if sig_pref == "all":
                    write_overlap(chemicals_file="%s/chemicals.txt" % (INPUTSPREFIX), sig_pref=sig_pref, forced=opts.forcesig)
                else:
                    # write both the sig and non-sig chemicals overlap
                    write_overlap(chemicals_file="%s/sig-chemicals.txt" % (INPUTSPREFIX), sig_pref=sig_pref, forced=opts.forcesig)
                    nonsig_pref = "nonsig-%s-%s-%s" % (scope, opts.sig_cutoff_type, str(opts.sig_cutoff).replace('.','_'))
                    write_overlap(chemicals_file="%s/nonsig-chemicals.txt" % (INPUTSPREFIX), sig_pref=nonsig_pref, forced=opts.forcesig)

            if opts.overlap_graph:
                # this also uses the sig-chemicals file, I just let the function handle it
                post_overlap_graph_to_graphspace(
                    interactome, opts.overlap_graph, total_num_chemicals,
                    num_sig_chemicals=len(curr_chemicals), sig_pref=sig_pref, forced=opts.forcesig
                )
# --------------------------------------------------
# -------------------- END MAIN --------------------


def generateResponseNetworks(chemicals, interactome, max_k=1000, pathlinker=False, k=500,
                             edge_penalty=None, rec_tfs_penalty=False, forced=False):
    """ Generate a response network for each chemical by running EdgeLinker (or PathLinker) using the set of responsive receptors and TFs as the sources and targets 
    *chemicals*: the chemicals to generate response networks for
    *max_k*: the total number of paths to write for edgelinker output
    *pathlinker*: run pathlinker instead of edgelinker
    *k*: number of paths to output for pathlinker
    *forced*: overwrite the output files if they already exist
    """
    # run pathlinker on every chemical, 
    print("Creating response networks for %d chemicals using the interactome '%s'" % (len(chemicals), interactome))
    # TODO implement the group run so reading the interactome (bottleneck) is only done once
    #t_utils.checkDir(out_dir)
    #run_edgelinker.runGroupEdgeLinker(interactome, rec_tfs_file, "%s/%s"%(out_dir, chemical), max_k)
    if pathlinker:
        print("Using PathLinker to generate the response networks")
    else:
        print("Using EdgeLinker to generate the response networks")

        if edge_penalty is not None:
            print("\tedge penalty: %s" % (str(edge_penalty)))
        if rec_tfs_penalty is True:
            print("\tapplying penalty to sources and targets (receptors and tfs)")

        out_dir = "%s/edgelinker"%(RESULTSPREFIX)
        t_utils.checkDir(out_dir)
        out_pref = "%s/" % out_dir
        edgelinker_in_files, edgelinker_out_files = run_edgelinker.createGroupEdgeLinker(chemicals, INPUTSPREFIX, out_pref)
        print("Running edgelinker on %d chemicals: %s, %s" % (len(chemicals), edgelinker_in_files, edgelinker_out_files))
        run_edgelinker.runEdgeLinker(interactome, edgelinker_in_files, edgelinker_out_files, max_k,
                     edge_penalty=edge_penalty, rec_tfs_penalty=rec_tfs_penalty, multi_run=True)

        # check to see if any of them failed, and output the failed runs
        for chemical in chemicals:
            results_file = "%s/%s-paths.txt"%(out_dir, chemical)
            if not os.path.isfile(results_file) or os.path.getsize(results_file) < 50:
                sys.stderr.write("ERROR: failed to generate response network for chemical: %s\n" % (chemical))
                # cleanup the errored out files
                results_file2 = "%s/%s-ranked-edges.txt"%(out_dir, chemical)
                os.remove(results_file)
                os.remove(results_file2)

        # the rest is for running pathlinker so stop here
        return 

    for chemical in tqdm(chemicals):
        # if version is netpath or kegg, use the different type of input file
        if VERSION in PATHLINKERDATAVERSIONS: 
            interactome = INTERACTOMES[VERSION] % chemical 
            #print("using interactome %s" % interactome)

        #chemical = chemicals[i]
        #print("\n----- %s %s %d -----"%(CHEMICAL_MAP.chemIDtoName[chemical], chemical, i))
        #interactomes = [opts.interactome]
#            chemical_weighted_interactomes = []

        #if opts.change_weight:
        #    pass
#                for change_weight in opts.change_weight:
#                #for change_weight in [0.3,0.1]:
#                    increase_weight = 1 + change_weight
#                    decrease_weight = 1 - change_weight
            # TODO change the dir to be in the inputs
#                    #ids_file = "results/chemicals/change_weights/%s-intermediate-acc.txt"%chemical
#                    t_utils.checkDir("results/chemicals/chemical_weights")
#                    interactome_name = opts.interactome.split('/')[-1]
#                    altered_interactome = "results/chemicals/chemical_weights/%s-%s-%s"%(chemical, str(increase_weight), interactome_name)
#                    write_altered_interactome(CHEMICAL_MAP, chemical, opts.interactome, altered_interactome, increase_weight, decrease_weight)
#                    chemical_weighted_interactomes.append((altered_interactome, change_weight))
        #else:
#            chemical_weighted_interactomes = [(interactome, '')]
#
#            # loop through the change weight options
#            for interactome, change_weight in chemical_weighted_interactomes:
#
#                if change_weight != '':
#                    # setup the output dir to be something like "weighted-1_10"
#                    change_weight = str(change_weight)[-1] + '_10'
#
#                #rec_tfs_file = "results/chemicals/pathlinker-inputs/%s-rec-tfs.txt"%chemical
        rec_tfs_file = REC_TFS_FILE % (INPUTSPREFIX, chemical)
        if pathlinker:

            out_dir = "%s/pathlinker"%(RESULTSPREFIX)
            t_utils.checkDir(out_dir)
            pathlinker_output_file = "%s/%sk_%s-ranked-edges.txt"%(out_dir, chemical, k)
            if not os.path.isfile(pathlinker_output_file) or forced:
                #runPathlinker(opts.interactome, rec_tfs_file, opts.k, out_dir+'/'+chemical, log_file='results/chemicals/log_files/%s.log'%chemical, ids=ids_file, cost_change=opts.change_weight)
                #runPathlinker(interactome, rec_tfs_file, opts.k, out_dir+'/'+chemical, log_file='', ids=ids_file, cost_change=opts.change_weight)
                runPathlinker(interactome, rec_tfs_file, k, out_dir+'/'+chemical, log_file='')

        #if opts.edgelinker:
        else:

            out_dir = "%s/edgelinker"%(RESULTSPREFIX)
            t_utils.checkDir(out_dir)
             
            edgelinker_output_file = "%s/%s-ranked-edges.txt"%(out_dir, chemical)
            edgelinker_paths_file = "%s/%s-paths.txt"%(out_dir, chemical)
            if not os.path.isfile(edgelinker_output_file) \
                    or os.path.getsize(edgelinker_paths_file) < 50 \
                    or forced:

                run_edgelinker.runEdgeLinker(interactome, rec_tfs_file, "%s/%s"%(out_dir, chemical), max_k)


def compute_stat_sig(scope, num_random_networks, k_to_test, forced=False):
    # after the random networks have been generated, run the statistical significance scripts to get p-values for each of the chemicals
    # example call:
    # python src/compute_stat_sig.py --compute-score-pval -n 1 1000 -k 10 -k 25 -k 50 -k 75 -k 100 -k 150 -k 200 -k 500 --out-dir outputs/2017_01-toxcast-signal-family-w0_15-p1_5-c1-no-bioseek/weighted/stats/stat-sig-permute-dir-undir/
    cmd = "python src/stat_sig/compute_stat_sig.py " + \
          " --chemicals %s/chemicals.txt " % (INPUTSPREFIX) + \
          " --compute-score-pval " + \
          " --num-random %d %d " % (num_random_networks[0], num_random_networks[1]) + \
          " -k %s " % (' -k '.join([str(k) for k in k_to_test])) + \
          " --out-dir %s/stats/stat-sig-%s/" % (RESULTSPREFIX, scope)
    if forced:
        cmd += " --forced "
    # TODO add parallelization to compute_stat_sig
    t_utils.runCommand(cmd)


def gen_figs(scope, k_to_test, forced=False):
    # TODO check to make sure the necessary files (such as the statistical significance) have been generated first
    commands_to_run = []

    # create some plots comparing the p-values to some other statistics
    cmd = "python src/plot/network/edge-weight-and-cost-hist.py " + \
          " --version %s " % (VERSION) + \
          " --compare-versions "
    commands_to_run.append(cmd) 

    cmd = "python src/plot/version_plots/network_summary_plots.py " + \
          " --version %s " % (VERSION) + \
          " --compare-versions "
    commands_to_run.append(cmd) 

#    cmd = "python src/plot/version_plots/k-type-hist.py " + \
#          " --version %s " % (VERSION) + \
#          " --compare-versions " + \
#          " -k %s " % (' -k '.join([str(k) for k in k_to_test])) 
#    commands_to_run.append(cmd) 
#
#    cmd = "python src/plot/version_plots/sig_plots.py " + \
#          "  --version %s " % (VERSION) + \
#          "  --compare-versions " + \
#          "  --scope %s " % (scope)
#    commands_to_run.append(cmd) 
#
#    # TODO fix this script
#    cmd = "python src/plot/version_plots/plot_degree.py " + \
#          " --version %s " % (VERSION) + \
#          " --compare-versions " + \
#          " --scope %s " % (scope)
#    #commands_to_run.append(cmd) 

    for cmd in commands_to_run:
        if forced:
            cmd += " --forced "
        print("\n")
        print("-"*40)
        t_utils.runCommand(cmd)


def write_overlap(chemicals_file, sig_pref="all", forced=False):
    # write the current chemical set to a file
    command = "python src/plot/overlap/plot_overlap_heatmap.py " + \
              " --chemicals %s " % (chemicals_file) + \
              " --out-dir  %s/stats/overlap/%s/ " % (RESULTSPREFIX, sig_pref) + \
              " --paths-dir %s/edgelinker/ " % (RESULTSPREFIX) + \
              " --write-overlap " + \
              " --plot-overlap " + \
              " --plot-hits-overlap "
              #" --inputs-dir %s " % (INPUTSPREFIX)
#              " --overlap-pvals "

    out_file = "%s/stats/overlap/%s/overlapping-proteins.txt" % (RESULTSPREFIX, sig_pref)
    if forced is True or not os.path.isfile(out_file):
        t_utils.runCommand(command)
    else:
        print("\t%s already exists. Use --forcesig to overwrite it" % (out_file))

    # also run it with the family nodes split
    command += " --split-family-nodes "
    out_file = "%s/stats/overlap/%s/split-overlapping-proteins.txt" % (RESULTSPREFIX, sig_pref)
    if forced is True or not os.path.isfile(out_file):
        t_utils.runCommand(command)
    else:
        print("\t%s already exists. Use --forcesig to overwrite it" % (out_file))


def post_overlap_graph_to_graphspace(interactome, overlap_cutoff, num_chemicals, num_sig_chemicals=None, sig_pref="all", forced=False):
    """
    Post the graph containing the nodes and edges common to > overlap_cutoff of the networks
    *chemicals_file*:
    *overlap_cutoff*: Number between 0 and 1. Generally 0.5
    *sig_pref*: scope-sig_cutoff_type-sig_cutoff. Used to get the overlap of the significant networks. Default "all" is the overlap of all networks
    *num_sig*: number of significant chemicals
    """

    if sig_pref == "all":
        chemicals_file = "%s/chemicals.txt" % (INPUTSPREFIX)
        graph_name = "overlap-graph-%s-0_5-color-%s" % (num_chemicals, VERSION)
        params = [(chemicals_file, sig_pref, graph_name)]
    else:
        # also post the non-significant overlap graph 
        sig_graph_name = "overlap-graph-sig%s-0_5-color-%s" % (num_sig_chemicals, VERSION)
        nonsig_graph_name = "overlap-graph-nonsig%s-0_5-color-%s" % (num_chemicals-num_sig_chemicals, VERSION)
        params = [("%s/sig-chemicals.txt"%(INPUTSPREFIX), sig_pref, sig_graph_name),
                  ("%s/nonsig-chemicals.txt"%(INPUTSPREFIX), "nonsig-"+sig_pref, nonsig_graph_name)]
    ev_version = t_utils.get_ev_version(VERSION)

    # TODO add this file to toxcast_settings.py(?)
    ctd_support_file = "/data/jeff-law/data/CTD/CTD_chem_gene_ixns.tsv"

    for chemicals_file, sig_pref, graph_name in params:
        node_overlap_file = "%s/stats/overlap/%s/overlapping-proteins.txt" % (RESULTSPREFIX, sig_pref)
        edge_overlap_file = "%s/stats/overlap/%s/overlapping-edges.txt" % (RESULTSPREFIX, sig_pref)
        # TODO update the evidence version
        command = "python src/graphspace/post_to_graphspace_overlap.py " + \
                " --ppi %s " % (interactome) + \
                " --version %s" % (ev_version)+ \
                " --sourcetarget %s/rec-tfs/all-rec-tfs.txt" % (INPUTSPREFIX) + \
                " --edge-overlap %s" % (edge_overlap_file) + \
                " --node-overlap %s" % (node_overlap_file) + \
                " --sig-chemicals %s" % (chemicals_file) + \
                " --overlap-cutoff %s" % (str(overlap_cutoff)) + \
                " --graph-name %s" % (graph_name) + \
                " --tag %s --tag overlap-graph --tag %s" % (VERSION, sig_pref) + \
                " --ctd-support %s" % (ctd_support_file)
                # TODO
                # "--group ToxCast2 " + \
        if VERSION in t_settings.EVIDENCE_FILES:
            evidence_file = t_settings.EVIDENCE_FILES[VERSION]
            command += " --ev-file %s" % (evidence_file)

        # This file is written by the script. If it already exists, then skip posting to graphspace
        ctd_support_amount = "%s-ctd-support.txt" % (node_overlap_file)
        if forced is True or not os.path.isfile(ctd_support_amount):
            t_utils.runCommand(command)
        else:
            print("\t%s already exists. Use --forcesig to post the overlap graph to graphspace anyway" % (ctd_support_amount))


def write_altered_interactome(chemical, interactome, altered_interactome, increase_weight, decrease_weight):
    out = open(altered_interactome, 'w')
    hits = []
    non_hits = []
    for acc in CHEMICAL_MAP.chemicals[chemical]['accs']:
        if acc not in CHEMICAL_MAP.chemicals[chemical]['rec'] and acc not in CHEMICAL_MAP.chemicals[chemical]['tfs']:
            if CHEMICAL_MAP.chemicals[chemical]['accs'][acc]['hit'] > 0:
                hits.append(acc)
            elif CHEMICAL_MAP.chemicals[chemical]['accs'][acc]['hit'] == 0:
                non_hits.append(acc)

    with open(interactome, 'r') as edges:

        for line in edges:
            if line[0] == '#':
                out.write(line)
            else:
                edge = line.rstrip().split('\t')
                u = edge[0]
                v = edge[1]
                w = float(edge[2])
                evidence = edge[3]

                if u in non_hits or v in non_hits:
                    w = w * decrease_weight

                elif u in hits or v in hits:
                    w = w * increase_weight
                    if w > 0.9:
                        w = 0.9

                out_line = '\t'.join([u,v,str(w),evidence]) + '\n'
                out.write(out_line)
    out.close()


def runPathlinker(interactome, rec_tfs_file, k, output_prefix, log_file='', ids='', cost_change=''):
    #sys.stdout.write("Runnin pathlinker on %s\n"%rec_tfs_file)
    sys.stdout.write("Runnin pathlinker on %s\n"%output_prefix)
    command = "python /home/jeffl/git-workspace/PathLinker/PathLinker.py %s %s -k %s -o %s --allow-mult-targets --allow-mult-sources --write-paths "%(interactome, rec_tfs_file, k, output_prefix)
    if ids and cost_change:
        command += "--tested-assays %s --increase-weight %0.1f --decrease-weight %0.1f"%(ids, 1+cost_change, 1-cost_change)
    if log_file:
        command += "> %s 2>&1"%(log_file)
     
    t_utils.runCommand(command, error_message="ERROR: Pathlinker failed on %s!\n"%rec_tfs_file)


def edgelinkerPermuteNetworks(interactome, scope, num_random_networks, max_k, k_to_test, forcealg, 
                             parallel_processes, super_computer, num_iterations=10,
                             printonly=False, cleanup=False, single_chem=[], send_email=None):
    """
    *num_iterations*: # edge swaps to perform (num_iterations * # edges in interactome)
    """

    out_dir = "%s/random-%s" %(RESULTSPREFIX, scope)
    t_utils.checkDir(out_dir)
    orig_interactome = interactome

    # just use the full interactome for now
    if VERSION in ["netpath-pathlinker-signaling-children-reg", "kegg-pathlinker-signaling-children-reg"]:
        # TODO use each pathway's interactome?
        interactome = INTERACTOMES["2015pathlinker"]

    # directory and file template for the stats of the generated random networks 
    stats_dir = "%s/stats/stat-sig-%s" % (RESULTSPREFIX, scope)
    rand_scores_k = "%s/rand-networks/rand-%%d-med-scores-k.txt" % (stats_dir)

    # check to make sure the chemical-k-median-scores.txt file used for summarizing the random networks exists
    # UPDATE 2017-01-04: use the median path score of the k+ties paths rather than the score of the kth path
    chemical_k_scores = "%s/chemical-k-median-scores.txt" % (stats_dir)
    # modify the paths a bit if this is running on baobab
    if super_computer:
        # TODO if the interactome has already been copied, don't need to copy it again
        # for now, just use the same chemical_k_scores
        if forcealg or not os.path.isfile(chemical_k_scores):
            # Copy the interactome to the baobab compute nodes
            t_utils.copyToBaobabNodes(interactome)
        else:
            print("Assuming interactome already copied to baobab compute nodes: %s\n\tIf not, use --forcealg to ensure it gets copied" % (interactome))
        t_utils.checkDir("%s/qsub_files" % out_dir)
        # keep the log files on /data
        qsub_out_dir = os.path.abspath(out_dir)
        # if cleanup is specified, write the results to the baobab nodes, and then delete them there
        # otherwise the results will be written to their original location
        if cleanup:
            # TODO this is a temporary fix for writing output to baobab nodes
            # 2016-11-15: move the file-writing and reading to the /localdisk partition on the baobab nodes (baobab-1, baobab-2, etc)
            # all of the reading and writing is slowing down /data (holmes) a lot
            out_dir = re.sub("^/data", "/localdisk", os.path.abspath(out_dir))
            interactome = re.sub("^/data", "/localdisk", os.path.abspath(interactome))
        else:
            print("\nTip: Use the --cleanup to write the files to individual baobab nodes rather than /data (holmes) offering a big speedup")

    # check to make sure the chemical-k-median-scores.txt file used for summarizing the random networks exists
    if forcealg or not os.path.isfile(chemical_k_scores):
        # Run the command to write the file if it does not exist
        # Whatever k are specified here will be used for all of the random networks (by parsing the header line of the file)
        # the user should select which k the significance values will be computed at
        # I am now using the seven values -k 25 -k 50 -k 75 -k 100 -k 150 -k 200 -k 500
        print("%s does not exist. Running compute_stat_sig.py with the --write-counts option to write it." % (chemical_k_scores))
        cmd = "python src/compute_stat_sig.py " + \
              " --chemicals %s/chemicals.txt " % (INPUTSPREFIX) + \
              " --k-limit %d " % (max_k) + \
              " --group-by-prob " + \
              " --write-counts " + \
              " --paths-dir %s/edgelinker " % (RESULTSPREFIX) + \
              " --out-dir %s " % (stats_dir) + \
              " -k %s " % (' -k '.join([str(k) for k in k_to_test]))
              #" -k 20 -k 50 -k 200 -k 500 "
        t_utils.runCommand(cmd)

    # command for permute
    command_template = "time python src/permute/permute_and_edgelinker.py \
            --version %s \
            --interactome %s \
            --inputs-dir %s \
            --out-dir %s --k %d \
            --num-iterations %d "\
            % (VERSION, interactome, INPUTSPREFIX, out_dir, max_k, num_iterations)
    command_template += " --write-score-counts %s " % (stats_dir)
    if scope == 'permute-dir':
        # if no options are specified, swapping will be all directed
        pass
    if scope == 'permute-undir':
        command_template += " --undirected "
    if scope == 'permute-dir-undir' or scope == 'permute-dir-undir-wb':
        command_template += " --swap-phys-sig-sep "
    if scope == 'permute-dir-undir-wb':
        command_template += " --split-by-weight 0.6 "
    if scope == 'local':
        all_rec_tfs_file = "%s/rec-tfs/all-rec-tfs.txt" % (INPUTSPREFIX)
        command_template += " --permute-rec-tfs %s " % (all_rec_tfs_file)
    if scope == 'global':
        all_rec_tfs_file = "%s/random-%s/all-rec-tfs.txt" % (INPUTSPREFIX, scope)
        global CHEMICAL_MAP
        # load the chemicals map
        if CHEMICAL_MAP is None: 
            CHEMICAL_MAP = t_utils.loadChemicalMap(orig_interactome)
        # write the set of all human rec and tfs
        CHEMICAL_MAP.write_global_rec_tfs(
            "%s/inputs/pathlinker-data/human-rec-tfs.txt" % (os.getcwd()),
            all_rec_tfs_file)
        command_template += " --permute-rec-tfs %s " % (all_rec_tfs_file)
    if forcealg:
        command_template += " --forced "
    if cleanup:
        command_template += " --cleanup "
    if single_chem is not None:
        command_template += " ".join([" --single-chem %s" % (chem) for chem in single_chem])

    # split up the jobs to run by the number of parallel processes
    # for example, running 1 to 100 jobs with 24 cores would be range(1,101,24) which returns [1, 25, 49, 73, 97]
    # 2016-11-17 Rather than loop through every set of indexes, just check if the desired output file exists or not
    jobs_to_run = []
    for i in range(num_random_networks[0], num_random_networks[1]+1):
        # if the stats for this index do not already exist, add this index to the job list
        if forcealg or not os.path.isfile(rand_scores_k % i):
            jobs_to_run.append(i)

    # parallel_processes*jobs_in_one_qsub jobs will be put in a single qsub script
    # these if else statements attempt to keep the number of jobs submitted to baobab to a minimum 
    # while also effectively using th eresources
    if len(jobs_to_run) > 3000:
        # there are only 6 baobab nodes, so if we had less than 600 jobs, 
        # we would not be using the resources effectively
        jobs_in_one_qsub = 50
    elif len(jobs_to_run) > 600:
        # if parralel_processes is 10, this put 10 permute_and_edgelinker.py jobs each with 10 random networks to make in one qsub, 
        # effectively putting 100 runs in one qsub script
        jobs_in_one_qsub = 10
    elif len(jobs_to_run) > 300:
        # up to 50 runs in one qsub script
        jobs_in_one_qsub = 5
    else:
        # up to 10 runs in 1 qsub script
        jobs_in_one_qsub = 1
    # first source the virtual environment
    # I'm not sure how to get this path automatically, so I'll just hard code it for now
    virtual_env = "source /data/jeff-law/tools/anaconda3/bin/activate  toxcast"
    # maybe I could just copy the entire current path(?)
    jobs = [virtual_env]
    start_index = 0
    for i in range(0, len(jobs_to_run), parallel_processes):
        # add the indexes of the jobs to run to the template
        end_index = i+parallel_processes if i+parallel_processes < len(jobs_to_run) else len(jobs_to_run)
        jobs.append(command_template + ''.join([" -r %d"%jobs_to_run[index] for index in range(i, end_index)]))

        if send_email and i+parallel_processes >= len(jobs_to_run):
            # if this is the last job, send an email indicating everything is finished
            # TODO also send an email if something fails?
            jobs.append(t_utils.setupEmailCommand(send_email, subject="ToxCast Permute Networks",
                                                  message="Finished `date`! Version: %s, scope: %s" % (VERSION, scope),
                                                  baobab=super_computer))

        # trick to group more jobs in the same qsub script
        # if this is not the last job and we have less jobs in the current command than we want, keep adding on to the command
        if len(jobs) < jobs_in_one_qsub and end_index < len(jobs_to_run):
            pass
        else:
            # skip running if printonly is specified
            if printonly:
                print(jobs)
                continue
            if super_computer:
                # run it on the super computer!
                qsub_pref = "%s/qsub_files/%d-%d" % (qsub_out_dir, jobs_to_run[start_index], jobs_to_run[end_index-1])
                #submit_jobs_super_computer([command], qsub_pref, ppn=parallel_processes)
                # 2016-12-01: Right now I'm running 10 jobs on a single node with 24 cores because
                # running more than 12 at the same time started causing jobs to quit with a -9 exit code (memory limit exceeded)
                submit_jobs_super_computer(jobs, qsub_pref, ppn=24)
            else:
                # run the command locally. It will start child processes equal to the specified number of parallel processes
                # TODO add local parallelization
                for job in jobs:
                    t_utils.runCommand(job)
                    #subprocess.check_call(job.split())

            # start the list of jobs over
            jobs = [virtual_env]
            # the new start index is the current end_index
            start_index = end_index


def submit_jobs_super_computer(jobs, qsub_pref, job_template='', nodes=1, ppn=24, walltime='15:00:00'):
    """ Function to write and submit a qsub bash script to the baobab cluster.
    Should only be run on the baobab cluster.
    This wiki page has the info for getting setup on baobab: http://wiki.cs.vt.edu/wiki/Cbb/Baobab
    *jobs* a list or set of commands/jobs that would normally be called with subprocess.check_output
    *qsub_pref* prefix to qsub file and logs (-out.txt and -err.txt)
    *job_template* a template of a job that has bash variables in it to be written to a loop. jobs should then be a list of (string) numbers or something
    *nodes* = total number of nodes you need
    *ppn* = processors per node that you will need
    *walltime* = amount of time your job will be allowed before being forcefully removed. 'HH:MM:SS'
    """
    qsub_file = "%s-qsub.sh" % qsub_pref
    name = qsub_pref.split('/')[-1] + '-' + VERSION
    # start the qsub file 
    with open(qsub_file, 'w') as out:
        out.write('#PBS -l nodes=%d:ppn=%d,walltime=%s\n' % (nodes, ppn, walltime))
        # set the job name
        out.write('#PBS -N %s\n' % (name))
        out.write('#PBS -o %s-out.txt\n' % (qsub_pref))
        out.write('#PBS -e %s-err.txt\n' % (qsub_pref))
        out.write('#####################################################################################\n')
        # edgelinker only runs in the directory the src code is in for some reason
        out.write('cd %s\n' % (os.getcwd()))
        # make sure the class path is set for the java command line parser
        out.write('export CLASSPATH="/home/jeffl/git-workspace/EdgeLinker-fork/commons-cli-1.3.1.jar:$CLASSPATH"\n')
        out.write('echo "Job Started at: `date`"\n')
        if not job_template:
            # write the jobs, as well as an echo (print) in a bash for loop to the qsub file
            out.write('\n'.join(['echo "%s"\n%s' % (cmd, cmd) for cmd in jobs]) + '\n')
        else:
            for_loop = 'COUNTER="%s"\n' % (' '.join(jobs))
            for_loop += 'for i in $COUNTER; do\n'
            for_loop += '\techo %s\n' % (job_template)
            for_loop += '\t%s\n' % (job_template)
            for_loop += 'done\n'
            out.write(for_loop)
        out.write('echo "Job Ended at: `date`"\n')
        # TODO some kind of email or notification if any of the jobs failed
    #command = "qsub -v interactome=%s,rec_tfs_file=%s,output_prefix=%s,max_k=%s src/run_edgelinker_baobab.sh" % (interactome, rec_tfs_file, out_dir+"/%d-random"%i, max_k)
    cmd = "qsub %s" % (qsub_file) 
    # run the qsub command from the system
    # if this is on baobab, then run as normal
    if 'baobab' in socket.gethostname():
        t_utils.runCommand(cmd)
    else:
        # otherwise run over ssh
        # have to give the path to the executable because the PATH variable isn't loaded when running over ssh
        t_utils.runCommandOnBaobab(cmd) 


def runCommandParallel(processes, max_processes, command, hide_output=False):
    # TODO Popen doesn't let you check the success/failure of the command. 
    # How do you know if it was successful or not?
    if hide_output:
        processes.add(subprocess.Popen(command.split(), stdout=open('/dev/null', 'w'), stderr=open('/dev/null', 'w')))
    else:
        processes.add(subprocess.Popen(command.split()))
    #print(len(processes), max_processes)
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update([
            p for p in processes if p.poll() is not None])
    return processes


def parseArguments(args):
    usage = 'master-script.py [options]\n'
    parser = OptionParser(usage=usage)

    ## Common options
    parser.add_option('','--version',type='string', action='append',
            help="Version of the PPI to run. Can specify multiple versions and they will run one after the other. Options are: %s." % (', '.join(t_settings.ALLOWEDVERSIONS)))
    parser.add_option('','--response-networks',action='store_true', default=False,
            help='Generate the response networks by running edgelinker with the set of perturbed receptors and tfs for each chemical. Use --pathlinker to run pathlinker instead')   
    parser.add_option('-S','--single-run', type='string', action='append',
            help='Run only a single chemical. Can specify multiple chemicals by listing this option multiple times.\n\tUse "file" to read the list of chemicals from the inputs/version/single_chem.txt file')   
    parser.add_option('','--forceinputs',action='store_true', default=False,
            help='Force overwriting the parsed input files')   
    group = OptionGroup(parser,'Algorithms')
    group.add_option('-P','--pathlinker',action='store_true', default=False,
            help='Run PathLinker instead of edgelinker on the set of chemicals using perturbed receptors and tfs')   
    #group.add_option('','--k',type='int',default='200',
    #        help='k to use when running pathlinker\n\tDefault = 200')   
    group.add_option('','--edgelinker-k',type='int',default='1000',
            help='maximum number of paths to keep when running edgelinker. Useful to not generate too much data (~2TB if all paths of each random network are kept)\n\tDefault = 1000')   
    group.add_option('-f','--forcealg',action='store_true', default=False,
            help='Force the algorithm to run even if the output is there already. Also works for forcing random network generation')   
    parser.add_option_group(group)

    group = OptionGroup(parser,'GraphSpace')
    group.add_option('-g','--graphspace',action='store_true', default=False,
            help='post response networks to graphspace')   
    group.add_option('','--k-to-post', type='int', default=200,
            help='Value of k to test for significance. Multiple k values can be given.')
    group.add_option('','--forcepost',action='store_true', default=False,
            help='Force the network to be posted to graphspace even if json file already exists')   
    group.add_option('-p','--postfix',action='store',
            help='add the postfix to the end of the chemical name to compare different options')   
    group.add_option('', '--revigo-colors', action='store_true', default=False,
            help="add colors to the nodes according to which REVIGO function/goterm it belongs to. GO terms must be manually selected based on REVIGO results")
    parser.add_option_group(group)

    group = OptionGroup(parser, 'Statistics')
    group.add_option('','--stat-sig', action='store_true',
            help='Compute the statistical significance for given k values (specified by --k-to-test).\n\tUses the --scope and --num-random-networks options')
    group.add_option('','--sig-cutoff', type='float',
            help='Significance (p-value) cutoff for selecting chemicals. Uses the --scope option.')
    group.add_option('','--sig-cutoff-type', type='string', default='FDR',
            help='Type of significance value to use. Default="FDR".\n\tValid options: "BF" (Bonferonni-corrected), "FDR" (False discovery rate q-value)')
    group.add_option('','--overlap', action='store_true', default=False,
            help='Write the # of times each protein is found in the response networks. Can be used with --sig-cutoff to get the overlap of the significant networks')
    group.add_option('','--overlap-graph', action='store', type="float", 
            help='Post the graph containing the nodes and edges common to > overlap_graph (cutoff) of the networks.' +
                 'Can be used with --sig-cutoff to get the overlap of the significant networks')
    group.add_option('','--gen-figs',action='store_true', default=False,
            help='Generate figures about various aspects of the networks including edge types, weights, statistical significance, and more')   
    group.add_option('','--forcesig',action='store_true', default=False,
            help='Force the scripts to compute statistical significance even if the output is there already')   
    parser.add_option_group(group)

    group = OptionGroup(parser,'Random Network Generation')
    #group.add_option('-R','--permute_network',action='store_true',
    #        help='Run edgelinker on a permuted network using the receptors and tfs for each chemical. Will be used to compute statistical significance of networks')   
    group.add_option('-r','--random-networks',action='store_true', default=False,
            help='Run pathlinker on random networks or random sets of receptors and tfs for each chemical. Will be used to compute statistical significance of networks')   
    group.add_option('','--scope',action='append', 
            help="type of permutation test to run. options are: 'local', 'global' and 'permute'." +
            "\n\t'local' is the set of assayed rec and tfs. \n\t'global' is the entire set of human rec and tfs. " +
            "\n\t'permute-undir' will run edgelinker on a permuted network where all edges are undirected using the receptors and tfs for each chemical." +
            "\n\t'permute-dir' will run edgelinker on a permuted network where all edges are directed using the receptors and tfs for each chemical." +
            "\n\t'permute-dir-undir' will permute the signaling and physical interactions (directed and undirected edges) separately." +
            "\n\t'permute-dir-undir-wb' will also permute the edges above and below a weight of 0.6 separately (wb = weighted bins).")   
    group.add_option('-k','--k-to-test', type='int', action='append',
            help='Value of k to test for significance. Multiple k values can be given. Suggested: 25, 50, 75, 100, 150, 200, 500')
    group.add_option('-n','--num-random-networks',type='int',default=(1,100),nargs=2,
            help="# of random networks to generate and run pathlinker on. Specify 2 integers (separated by a space): start_index end_index. \n\tDefault = 1 100")   
    #group.add_option('','--chemical-range',action='store', type='int', nargs=2,
    #        help='range of chemicals to run the random sets on. Requires two integers. For example: --chemical-range 100 200. Reverse list also supported (i.e. 200 100). Useful for generating random networks in parallel.')   
    group.add_option('','--parallel-processes',action='store', type='int', default=1,
            help='# of processes to run in the background simultaneously. Default = 1')   
    group.add_option('','--super-computer',action='store_true', default=False,
            help='Generate random networks using the baobab cluster. Can only be run on a lab computer or from baobab (ssh baobab.cbb.lan)')   
    group.add_option('','--cleanup',action='store_true', default=False,
            help='Delete edgelinker output the intermediate files generated. Also writes the files to individual baobab nodes rather than /data (holmes) offering a big speedup (holmes gets really slow with all of the I/O). \n\tCurrently only implemented for random network generation')
    # Was used for the 'local' and 'global' settings when random sets of receptors and tfs were chosen for each chemical
    group.add_option('','--group-num-rectfs',action='store_true', default=False,
            help='For the \'local\' and \'global\' scopes, group the random network generation by the # of rec and tfs ' + \
                    'rather than random networks for each chemical. Cuts data generation ~1/2 for 464 chemicals.')   
    group.add_option('','--send-email',action='store', type='string',
            help='Send an email (to specified recipient, from jeffreynlaw@gmail.com) after the random networks have finished generating')
    parser.add_option_group(group)

    group = OptionGroup(parser, 'Other Options')
    group.add_option('','--printonly',action='store_true', default=False,
            help='Just print the commands that would be run. Currently only implemented for random network generation')
    parser.add_option_group(group)

    # previous options
#    parser.add_option('-o','--out-dir',action='store', default="results/chemicals/pathlinker-outputs",
#            help='output dir to store pathlinker results for each chemical')   
#    parser.add_option('','--change-weight',action='append', type='float',
#            help='Change the weight of  edges out of and into acc by this amount')   
#    parser.add_option('','--overlap',action='store_true', \
#            help='compute the pval of overlap of chemical A to B')   
#    parser.add_option('','--compute-prec-rec',action='store_true', \
#            help='Compute the precision recall values for the chemicals')   
    #parser.add_option('','--write-paths',action='store_true',default=False, help='')

    # parse the command line arguments
    (opts, args) = parser.parse_args()

    #if not opts.chemicals and not opts.random_networks and not opts.graphspace and not opts.overlap:
    #    print("Must specify to run pathlinker on the set of chemicals or the random set of receptors")
    #    sys.exit(1)

    #if opts.chemicals and not (opts.pathlinker or opts.edgelinker):

    for version in opts.version:
        if version not in t_settings.ALLOWEDVERSIONS:
            print("ERROR: version '%s' not an allowed version. Options are: %s." % (version, ', '.join(t_settings.ALLOWEDVERSIONS)))
            sys.exit()

    if opts.scope is not None:
        for scope in opts.scope:
            if scope not in ['local', 'global', 'permute-dir', 'permute-dir-undir', 'permute-dir-undir-wb', 'permute-undir']:
                print("ERROR: %s not a valid option for scope. Valid options are: %s" %(scope, ', '.join(['local', 'global', 'permute-dir', 'permute-dir-undir', 'permute-dir-undir-wb', 'permute-undir'])))
                sys.exit(1)
    # This could be fine if the user wanted to generate them in reverse order (i.e. 464 -> 1) to try to eliminate problems of running edgelinker on the same chemicals
#    if opts.chemical_range:
#        if opts.chemical_range[0] > opts.chemical_range[1] or opts.chemical_range[0] < 0:
#            print("ERROR: chemical range is invalid. First number must be < second number. Given: %d, %d"%(opts.chemical_range[0],opts.chemical_range[1]))
#            sys.exit(1)

    return opts

if __name__ == '__main__':
    main(sys.argv)
