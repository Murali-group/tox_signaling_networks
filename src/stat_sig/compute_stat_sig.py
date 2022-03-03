#!/usr/bin/env python
""" Script for computing the statistical significance of each chemical's network compared to the generated random networks
"""

from optparse import OptionParser
import os, sys, traceback
from src import toxcast_utils as t_utils
from src.utils import file_utils as utils
# utilize the gpdPerm script to estimate the generalized Pareto distribution
from src.stat_sig.gpdPerm import gpdPerm
from tqdm import tqdm
import numpy as np
import gc


class StatSig:

    def __init__(self, paths_dir=None, out_dir=None, random_paths_dir=None, k_limit=200, num_random=(1,10000), group_by_prob=False):
        self.k_limit = k_limit
        self.num_random = num_random[1]
        self.num_random_sets = num_random
        self.out_dir = out_dir
        self.group_by_prob = group_by_prob
        self.pval_cutoff = 0.05

        # define the paths to the pathlinker output files
        if paths_dir is not None:
            self.pathlinker_paths = paths_dir+'/%s-paths.txt'
            self.pathlinker_ranked_edges = paths_dir+'/%s-ranked-edges.txt'
        if random_paths_dir is not None:
            self.random_pathlinker_paths = random_paths_dir+'/%s/%d-random-paths.txt'
            self.random_pathlinker_ranked_edges = random_paths_dir+'/%s/%d-random-ranked-edges.txt'

    def getProteinCount(self, ranked_edges_file, k_limit):
        """ Get the number of proteins or nodes in the top k_limit paths 
            *ranked_edges_file* is the file of ranked edges output by PathLinker
            *k_limit* is the path index at which to stop counting proteins
        """
        proteins = set()
        with open(ranked_edges_file, 'r') as file_handle:
            for line in file_handle:
                if line[0] == '#':
                    continue
                line = line.rstrip().split('\t')
                k = int(line[2])
                if k > k_limit:
                    break
                proteins.add(line[0])
                proteins.add(line[1])

        return len(proteins)

    def getPathStats(self, paths_file, k_limit, group_by_prob=False):
        """   Get the statistics of a single paths file 
            *paths_file* is the file output by CycLinker (output by Peter Steele's java implementation of cyclinker)
            *k_limit* int > 0 is the limit of the k # of paths to analyze
            returns a dictionary containing the path index as the key and the probability as the value
        """
        try:
            k_prob = {}
            #probs_count = {}
            with open(paths_file, 'r') as file_handle:
                for line in file_handle:
                    if line[0] != '#':
                        if len(line.split('\t')) >= 3:
                            # example pathlinker line:
                            #1   1.41162e-01 P00533|P42229|P10242
                            line = line.rstrip().split('\t')
                        else:
                            line = line.rstrip().split(' ')
                        k = int(line[0])
                        #path = line[2].split('|')
                        # limit the probability to 6 decimal places
                        probability = float("%0.06f" % float(line[1]))
                        if k > k_limit:
                            break

                        # count the number of paths that tie at this probability
                        #if group_by_prob:
                        #if probability not in probs_count:
                        #    probs_count[probability] = 0 
                        #probs_count[probability] += 1
                        #else:
                        k_prob[k] = probability
                            # move the probability back into the -log domain for adding the sums
                            #probabilities[k] = -1 * (math.log10(probability))
        except ValueError:
            print("\nFailed to get path stats for %s." % paths_file)
#            # If this is the final version, attempt to re-run cyclinker
#            # outputs/2016_05-pathlinker-signalling-5_100/weighted/random-global//2rec-9tfs/2633-random-paths.txt
#            if "2016_05-pathlinker-signalling-5_100" in paths_file and 'random-global' in paths_file:
#                print("Attempting to re-run cyclinker")
#                i = paths_file.split('/')[-1].split('-')[0]
#                out_dir = "%s/%s-random" % ('/'.join(paths_file.split('/')[:-1]), i)
#                interactome = "/home/jeffl/svnrepo/data/interactomes/human/2016_05/2016-06-02-human-ppi-mod-sig-5_100-weighted.txt"
#                rec_tfs_file = "inputs/2016_05-pathlinker-signalling-5_100/random-global/%s/%s-random-rec-tfs.txt" % (paths_file.split("/")[-2], i)
#                run_cyclinker.runCyclinker(interactome, rec_tfs_file, out_dir, k_limit + 1000)
#                # now try to get the pathStats again
#                # TODO could this create an endless loop?
#                k_prob, probs_count = self.getPathStats(paths_file, k_limit, group_by_prob)
#            else:
            raise
        if len(k_prob) == 0:
            print("\nFailed to get path stats for %s." % paths_file)
            #print("\nFailed to get path stats for %s. Quitting" % paths_file)
            #sys.exit()

        #return k_prob, probs_count
        return k_prob

    def computePathProbPval(self, chemical, chem_prob, rand_prob, k_limit, pval_sum=False, pval_min=False, reverse_pval=False):
        """   Compute the p-value of the chemical's path probabilities compared to the random network's probabilities
        *chem_prob* 
        *rand_prob* is a list of random network probabilities dictionaries
        *k_limit* is the limit of paths for which to sum and/or max
        *sum* is a T/F value to compute the sum of the 
        *min* is a T/F value to compute the minimum fraction of random paths with a probability >= the chemical's paths over any of the k-values
        """
        sum_prob = 0
        k_pvals = {}
        sum_pvals = {}
        sum_counts = {}
        rand_counts = {}
        #for k in range(1,min(len(chem_prob)+1, k_limit+1)):
        for k in sorted(chem_prob) if not reverse_pval else sorted(chem_prob, reverse=True):
            # troubleshoot missing k values
            #for i in xrange(len(rand_prob)):
            #    try:
            #        rand_prob[i][k]
            #    except KeyError:
            #        print("chemical %s is missing %d-random-paths.txt"%(chemical, i+1))
            #        sys.exit(1)
            # get the path probability at each k value, or the # of paths at each probability cutoff
            rand_prob_k = [rand_prob[i][k] if k in rand_prob[i] else 0 for i in xrange(len(rand_prob))]
            k_pvals[k] = t_utils.computePval(chem_prob[k], rand_prob_k, reverse_pval)

            if pval_sum:
                #pdb.set_trace()
                # compute the X_k pvalue
                sum_prob += chem_prob[k]
                rand_sum_probs = [] 
                for i in xrange(len(rand_prob)):
                    # The path probabilities of the random networks might not match up with the chemical's network
                    # Compute the sum of the # of paths above the probability cutoff (k) (for example: # of paths above a path probability of 0.5)
                    rand_sum_probs.append(sum([rand_prob[i][prob] for prob in rand_prob[i] if prob >= k]))
                sum_pvals[k] = t_utils.computePval(sum_prob, rand_sum_probs, reverse_pval)
                sum_counts[k] = sum_prob
                rand_counts[k] = rand_sum_probs

        return sum_counts, rand_counts, sum_pvals, k_pvals

    def plotCumulativeHisto(self, chem_prob, rand_prob, fig_file, k, pval, sum_k=False):
        """ plot the cumulative histogram of probabilities of the random network's k path
        with a dotted line showing the chemical's network k path probability
        """

        # here we want to use a reversed cumulative histogram
        fig, ax = plt.subplots()
        plt.hist(rand_prob, normed=1, histtype='step', cumulative=True)
        # us the regular cumulative for the -log space
        #plt.hist(rand_prob, normed=1, histtype='step', cumulative=-1)

        # plot the single chemical's probability line
        plt.axvline(chem_prob, linestyle='--')

        # set the labels and title 
        ax.set_ylabel('Fraction of Random Networks')
        if sum_k:
            ax.set_title("Cumulative Histogram of Sum Path Probability \nup to k%d pval=%0.3f"%(k, pval))
            ax.set_xlabel('sum of path probabilities')
        else:
            ax.set_title("Cumulative Histogram of Path Probability \nk%d pval=%0.3f"%(k, pval))
            ax.set_xlabel('path probability')

        plt.ylim(ymin=-0.05,ymax=1.05)

        # save the figure
        plt.savefig(fig_file)
        plt.close(fig)

    def get_chem_counts(self, k_list=None, use_counts=False):
        """ 
        *chem_scores_file*: 
            Default: "%s/chemical-k-scores.txt" % self.out_dir
        """
        if k_list is None:
            # use the default list of k values
            k_list = [10, 25, 50, 75, 100, 150, 200, 500]
        # the "use_counts" option doesn't have any effect because the counts column is the same in the two files.
        # Just the scores are different, but the scores aren't used with the --use-counts option
        if use_counts is True:
            chem_scores_file = "%s/chemical-k-scores.txt" % self.out_dir
        else:
            chem_scores_file = "%s/chemical-k-median-scores.txt" % self.out_dir
        # get the scores at each k in the already written file
        # mapping of a chemical to the k + ties counts for each k specified (in the file) 
        chem_counts = {}
        chem_scores = {}
        with open(chem_scores_file, 'r') as file_handle:
            header = file_handle.readline()
            for line in file_handle:
                if line[0] == '#':
                    continue
                line = line.rstrip().split('\t')
                chem = line[0]
                scores = [float(line[j]) for j in range(1, len(line), 2)]
                counts = [int(line[j]) for j in range(2, len(line), 2)]
                chem_scores[chem] = scores
                chem_counts[chem] = counts
                if len(counts) != len(k_list):
                    print("ERROR: # of k-to-test in chem-counts file does not match # of k-to-test specified  (%d != %d). Did you use the same --k-to-test (-k) for both?" % (len(counts), len(k_list)))
                    sys.exit(1)

        return chem_scores, chem_counts

    def get_rand_counts(self, chem_counts=None, use_counts=False, rand_out_template=None):
        """ Get the count of the # of paths at each score cutoff for each k from each random network scores file
        """
        if chem_counts is None:
            # the chem_counts dictionary is used to make sure the same # of k indices are used for both of the files
            chem_scores, chem_counts = self.get_chem_counts(use_counts=use_counts)
        if rand_out_template is None:
            if use_counts is True: 
                rand_out_template = "%s/rand-networks/rand-%%d-scores-k.txt" % (self.out_dir)
            else:
                rand_out_template = "%s/rand-networks/rand-%%d-med-scores-k.txt" % (self.out_dir)
        rand_counts = {}
        rand_scores = {}
        files_skipped = 0 
        # now get the random counts at each k 
        # mapping of a chemical to the dictionary of random network indices to the # of paths above the score cutoff from the chem_scores_file 
        for rand_index in tqdm(range(self.num_random_sets[0],self.num_random_sets[1]+1)):
            rand_counts_file = rand_out_template % rand_index
            if not os.path.isfile(rand_counts_file):
                # if the file was not generated by permute_and_cyclinker.py, then skip it
                tqdm.write("Warning: %s does not exist. Skipping it" % (rand_counts_file))
                files_skipped += 1
                continue
            with open(rand_counts_file, 'r') as file_handle:
                for line in file_handle:
                    if line[0] == '#':
                        continue
                    line = line.rstrip().split('\t')
                    chemical = line[0]
                    # UPDATE 2017-11-30: I am now tracking the scores as well. Keep the original functionality as they don't have the scores
                    if use_counts is True:
                        counts = [int(count) for count in line[1:]]
                    else:
                        scores = [float(line[j]) for j in range(1,len(line),2)]
                        counts = [int(line[j]) for j in range(2,len(line),2)]
                    if chemical not in rand_scores:
                        rand_scores[chemical] = {}
                    rand_scores[chemical][rand_index] = scores
                    if chemical not in rand_counts:
                        rand_counts[chemical] = {}
                    rand_counts[chemical][rand_index] = counts
                    if len(counts) != len(chem_counts[chemical]):
                        print("ERROR: # of k-to-test in random file != # of k-to-test in chem-counts file  (%d != %d). Did you use the same --k-to-test (-k) for both?" % (len(counts), len(chem_counts[chemical])))
                        print("\trand counts file: %s\n" % (rand_counts_file))
                        sys.exit(1)

        if files_skipped >= 1:
            print("\nA total of %d random counts files did not exist. Continuing\n" % (files_skipped))

        return rand_scores, rand_counts

    def compute_score_pval(self, k_list, use_counts=False, forced=False):
        k_pval_file = "%s/k-perm-pval.txt" % self.out_dir
        gpd_pval_file = "%s/gpd-pval.txt" % self.out_dir
        gpd_pvals_k200_file = "%s/gpd-pvals-k200.txt" % self.out_dir
        corr_p_qval = "%s/bfcorr_pval_qval.txt" % self.out_dir 
        if not forced and os.path.isfile(k_pval_file) and os.path.isfile(gpd_pval_file) and os.path.isfile(gpd_pvals_k200_file) and os.path.isfile(corr_p_qval):
            print("All output files (k-perm-pval.txt, gpd-pval.txt, gpd-pvals-k200.txt, bfcorr_pval_qval.txt) already exist. Use --forced to overwrite them.")
        else:
            chem_scores_file = "%s/chemical-k-scores.txt" % self.out_dir
            print("Getting the random counts above the scores of each chemical's network at k values in %s " % chem_scores_file)
            with open(chem_scores_file, 'r') as f: 
                header = f.readline().rstrip().split('\t')
            chem_scores, chem_counts = self.get_chem_counts(k_list, use_counts=use_counts)
            rand_scores, rand_counts = self.get_rand_counts(chem_counts, use_counts=use_counts)

            print("Computing the permutation and gpd p-values for each chemical")
            gpd_pvals = {}
            gpd_pvals_k200 = {} 
            perm_pvals = {}
            # compute the counts p-value  for each chemical
            for chemical in tqdm(chem_counts):
                perm_pvals[chemical] = {}
                gpd_pvals[chemical] = {}
                gpd_pvals_k200[chemical] = {}
                # and for each k
                # index of j in the list corresponds to the # of ties or the counts at each of the k values for every chemical
                for k_index in range(len(chem_counts[chemical])):
                    # comput ethe pvalue separately for each chemical and for each k specified 
                    if use_counts is False:
                        chem_count = chem_scores[chemical][k_index]
                        rand_chem_counts = [rand_scores[chemical][i][k_index] for i in rand_scores[chemical]]
                    else:
                        chem_count = chem_counts[chemical][k_index]
                        rand_chem_counts = [rand_counts[chemical][i][k_index] for i in rand_counts[chemical]]

                    perm_pvals[chemical][k_index] = t_utils.computePval(chem_count, rand_chem_counts, reverse=True)
                    # TODO How can I keep track of the goodness of fit p-value and such for each of the scores?
                    # compute the gpd pvalue. Will only be different for the permutation p-values less than 10 / # random networks (i.e. < 10/10000) 
                    gpd_pval, M, K, A, nexc, good_fit_pval = gpdPerm.est(chem_count, rand_chem_counts, Nexcmax=250)
                    gpd_pvals[chemical][k_index] = gpd_pval 
                    # keep track of all of the metrics for the k200.
                    # position 2 corresponds to k200 
                    # TODO I need a more permanent way to select k-200 
                    # example header: 
# #chemical   k-25-score  count   k-50-score  count   k-75-score  count   k-100-score count   k-150-score count   k-200-score count   k-500-score count
                    if '200' in header[k_index*2]:
                        gpd_pvals_k200[chemical][k_index] = (chem_count, gpd_pval, M, K, A, nexc, good_fit_pval)

            if not forced and os.path.isfile(k_pval_file):
                print("not writing %s. File already exists. Use --forced to overwrite it" % (k_pval_file))
            else:
                print("writing ", k_pval_file)
                # now write the output files
                with open(k_pval_file, 'w') as out:
                    out.write("#chemical\t" + '\t'.join([str(k) for k in k_list]) + '\n')
                    for chemical in sorted(perm_pvals):
                    #for k in k_list:
                        #out.write('\n'.join(["%s\t%d\t%0.4f\t%0.4f" % (chemical, chem_counts[chemical][perm_pvals[chemical][k] for chemical in perm_pvals])
                        out.write(chemical + '\t' + '\t'.join(["%0.4f" % (perm_pvals[chemical][j]) for j in perm_pvals[chemical]]) + '\n')

            if not forced and os.path.isfile(gpd_pval_file):
                print("not writing %s. File already exists. Use --forced to overwrite it" % (gpd_pval_file))
            else:
                print("writing ", gpd_pval_file)
                with open(gpd_pval_file, 'w') as out:
                    out.write("#chemical\t" + '\t'.join([str(k) for k in k_list]) + '\n')
                    for chemical in sorted(gpd_pvals):
                    #for k in k_list:
                        #out.write('\n'.join(["%s\t%d\t%0.4f\t%0.4f" % (chemical, chem_counts[chemical][perm_pvals[chemical][k] for chemical in perm_pvals])
                        out.write(chemical + '\t' + '\t'.join(["%0.4e" % (gpd_pvals[chemical][j]) for j in gpd_pvals[chemical]]) + '\n')

            # TODO write the q-value and/or the bonferroni-corrected p-value
            if not forced and os.path.isfile(gpd_pvals_k200_file):
                print("not writing %s. File already exists. Use --forced to overwrite it" % (gpd_pvals_k200_file))
            elif len(gpd_pvals_k200) == 0:
                print("not writing %s. No k of 200 found" % (gpd_pvals_k200_file))
            else:
                print("writing ", gpd_pvals_k200_file)
                with open(gpd_pvals_k200_file, 'w') as out:
                    out.write("#chemical\tk\tgpd pval\tM\tnexc\tk\ta\tgood_fit_pval\n")
                    for chemical in sorted(gpd_pvals_k200):
                        for j in gpd_pvals_k200[chemical]:
                            chem_count, gpd_pval, M, K, A, nexc, good_fit_pval = gpd_pvals_k200[chemical][j]
                            out.write("%s\t%d\t%0.4e\t%s\t%s\t%s\t%s\t%0.4f\n" % (chemical, chem_count, 
                                    gpd_pval, str(M), str(nexc), str(K), str(A), good_fit_pval))

            if not forced and os.path.isfile(corr_p_qval):
                print("not writing %s. File already exists. Use --forced to overwrite it" % (corr_p_qval))
            else:
                print("writing ", corr_p_qval)
                output = []
                num_chemicals = len(gpd_pvals)
                # Sort by the k200 BF corrected pvalue column
                for j in range(len(k_list)):
                    # get the chemical and the p-value in a tuple so they can be sorted
                    chem_pval = []
                    for chemical in sorted(gpd_pvals):
                        pval = gpd_pvals[chemical][j]
                        chem_pval.append((chemical, pval))
                    # sort the chemicals by the p-value
                    chem_corr_pvals = []
                    index = 1
                    for chemical, pval in sorted(chem_pval, key=lambda x: x[1]):
                        # compute the bonferroni corrected p-value by multiplying by the number of chemicals tested 
                        # limit the p-value at 1
                        bfcorr_pval = min(1, pval * num_chemicals)
                        # compute the FDR q-value by multiplying by the number of chemicals divided by the index
                        # also limit the q-value at 1
                        qval = min(1, pval * (num_chemicals / float(index)))
                        #print("pval, qval, num_chemicals, index", pval, qval, num_chemicals, index)
                        chem_corr_pvals.append((chemical, bfcorr_pval, qval)) 
                        index += 1
                    output.append(chem_corr_pvals)

                # now sort all of the lists by the order of the chemicals to get all of them on the same line for each chemical
                sorted_by_chemical_output = {}
                for chem_corr_pvals in output:
                    # sort by the first item which is the chemical
                    chem_corr_pvals = sorted(chem_corr_pvals)
                    for chemical, bfcorr_pval, qval in chem_corr_pvals:
                        if chemical not in sorted_by_chemical_output:
                            sorted_by_chemical_output[chemical] = []
                        sorted_by_chemical_output[chemical].append(bfcorr_pval)
                        sorted_by_chemical_output[chemical].append(qval)

                # sort the entire thing by the k200 bfcorr pval column (which is 6th column or column index 5)
                sorted_by_k200_bfcorr_pval = []
                for chemical in sorted_by_chemical_output:
                    sorted_by_k200_bfcorr_pval.append([chemical] + sorted_by_chemical_output[chemical])
                sorted_by_k200_bfcorr_pval = sorted(sorted_by_k200_bfcorr_pval, key=lambda x: x[5])

                with open(corr_p_qval, 'w') as out:
                    out.write("#chemical\t" + '\t'.join(["k%s-BFcorr-pval\tk%s-qval" % (k,k) for k in k_list]) + '\n')
                    for line in sorted_by_k200_bfcorr_pval:
                        pvals = '\t'.join(["%0.4e" % (p) for p in line[1:]])
                        out.write(line[0] + '\t' + pvals + '\n')
                    #out.write("\n".join([line[0] + '\t'.join(["%0.4e" pval for pval in line[1:]]) for line in sorted_by_k200_bfcorr_pval])

        # now print(some summary statistics)
        # TODO make this more autmoated
        perm_pvals = utils.readColumns(k_pval_file, *range(2,len(k_list)+2))
        corr_pvals = utils.readColumns(corr_p_qval, *range(2,len(k_list)*2+2,2))
        qvals = utils.readColumns(corr_p_qval, *range(3,len(k_list)*2+2,2))
        print("\n" + "-"*20)
        print("Computing the number of significant response networks at each k (corrected for the number of k tested)")
        #print("\talpha: %0.4f (%0.3f/%d)" % (self.pval_cutoff/len(k_list), self.pval_cutoff, len(k_list)))
        print("\talpha: %0.3f" % (self.pval_cutoff))
        print("pval_type\tsig_any_k\t%s" % ('\t'.join(["k%s"%k for k in k_list])))
        for pval_type, pvals in [('perm_pvals', perm_pvals), ('BF_gpd_pvals', corr_pvals), ('FDR_qvals', qvals)]:
            num_sig = [0]*len(k_list)
            sig_any_k = 0
            for chem_pvals in pvals:
                sig_chem = False
                for i in range(len(chem_pvals)):
                    # default of pval_cutoff is 0.05
                    #if float(chem_pvals[i]) < self.pval_cutoff/len(k_list):
                    if float(chem_pvals[i]) < self.pval_cutoff:
                        num_sig[i] += 1
                        sig_chem = True
                # if this chemical is significant at any of the k values, add 1 to the overall cont
                sig_any_k = sig_any_k + 1 if sig_chem else sig_any_k
            print("%s\t%d\t%s" % (pval_type, sig_any_k, '\t'.join([str(num) for num in num_sig])))

    def print_pval_summary(self, pval_lists, k_list):
        """ Function to print(out the # of networks with a p-value below the cutoff)
        *pval_lists* - a dictionary with a string of the type of p-value as the key and the list of # of networks passing the 
        *pval_list* - the list of 
        """ 

    def write_counts(self, chemicals_file, k_list, forced=False):
        """
        """
        print("Getting the score at each k, k+# of ties (count), and the median score of the k+ties paths for each k: " + str(k_list))

        self.chemicals = utils.readItemList(chemicals_file, col=1)
        t_utils.checkDir(self.out_dir)
        #chem_k_scores = {}
        
        # start with the header line
        out_string = "#chemical\t" + '\t'.join(['k-%d-score\tcount' % k for k in k_list]) + '\n'
        out_string_med = "#chemical\t" + '\t'.join(['k-%d-score-med\tcount' % k for k in k_list]) + '\n'

        #print("\t%s parsing chem values" % (chemical))
        for chemical in tqdm(sorted(self.chemicals)):
            # get the chemical's subnetwork path probabilities
            chem_paths = self.pathlinker_paths % chemical
            chem_k_score = self.getPathStats(chem_paths, self.k_limit, self.group_by_prob)

            #path_count = 0
            #sum_score_count = {}
            #for score in sorted(chem_score_count, reverse=True):
            #    path_count += chem_score_count[score]
            #    sum_score_count[score] = path_count

            # for each k in the list, get the corresponding path score
            #out_string += "\t".join([chemical] + ['%0.6f\t%d'%(chem_k_score[k],sum_score_count[chem_k_score[k]]) for k in k_list]) + '\n'
            # UPDATE 2017-12-15: for each k, get the median of the path scores up to that spot
            scores = list(chem_k_score.values())
            out_string_med += chemical
            out_string += chemical
            for k in k_list:
                ties_count = sum([1 for curr_k in chem_k_score if chem_k_score[curr_k] >= chem_k_score[k]])
                median_paths_score = np.median(scores[:ties_count])
                out_string_med += "\t%0.6f\t%d" % (median_paths_score, ties_count)
                out_string += "\t%0.6f\t%d" % (chem_k_score[k], ties_count)
            out_string += '\n' 
            out_string_med += '\n'

        # write both files so we have both the score at the given k, and the median score of the k+ties paths
        print("Writing %s/chemical-k-scores.txt" % (self.out_dir))
        with open("%s/chemical-k-scores.txt" % self.out_dir, 'w') as out:
            out.write(out_string)
        # I left the "count" column (k+# of ties) in this file just so they're formatted the same
        print("Writing %s/chemical-k-median-scores.txt" % (self.out_dir))
        with open("%s/chemical-k-median-scores.txt" % self.out_dir, 'w') as out:
            out.write(out_string_med)

    def write_rand_counts(self, chemicals=None, chemicals_file=None, k_list=None, forced=False):

        if k_list is None:
            # use the default list of k values
            k_list = [10, 25, 50, 75, 100, 150, 200, 500]

        if chemicals is None:
            self.chemicals = utils.readItemList(chemicals_file, col=1)
        else:
            self.chemicals = chemicals
        t_utils.checkDir("%s/rand-networks" % self.out_dir)
        rand_out_template = "%s/rand-networks/rand-%%d-med-scores-k.txt" % (self.out_dir)
       # chem_scores_file = "%s/chemical-k-scores.txt" % self.out_dir
#
#        print("Getting the random counts above the scores of each chemical's network at k values in %s and writing them to %s " % (chem_scores_file, rand_out_template))
#
#        # get the scores at each k in the already written file
#        chem_score_list = {}
#        with open(chem_scores_file, 'r') as file_handle:
#            for line in file_handle:
#                if line[0] == '#':
#                    continue
#                line = line.rstrip().split('\t')
#                score_list = []
#                for i in range(1, len(line), 2):
#                    score_list.append(float(line[i]))
#                chem_score_list[line[0]] = score_list
#                if len(score_list) != len(k_list):
#                    print("ERROR: # of k-to-test in chem-counts file does not match # of k-to-test specified  (%d != %d). Did you use the same --k-to-test (-k) for both?" % (len(score_list), len(k_list)))
#                    print("TODO: add the k-to-test option to permute_and_cyclinker.py")
#                    sys.exit(1)
        
        #print("\t%s parsing chem values" % (chemical))
        for rand_index in tqdm(range(self.num_random_sets[0],self.num_random_sets[1]+1)):
            # if the file already exists, then skip this index unless forced is specified
            if not forced and os.path.isfile(rand_out_template % rand_index):
                print("%s already exists. Skipping" % (rand_out_template % rand_index))
                continue
            #try:
            #out_string = "#chemical\t" + '\t'.join(['score-%0.6f' % score for score in score_list[chemical]]) + '\n'
            #out_string = ''
            out_string_med = ""
            for chemical in tqdm(sorted(self.chemicals)):
                # get the chemical's subnetwork path probabilities
                rand_paths = self.random_pathlinker_paths %(chemical, rand_index)
                rand_k_score = self.getPathStats(rand_paths, self.k_limit, self.group_by_prob)

#                rand_chem_count = {}
#                for score in chem_score_list[chemical]:
#                    # count the number of paths with a score greater than the chemical's score
#                    rand_chem_count[score] = sum([rand_score_count[score2] for score2 in rand_score_count if score2 >= score])
#
#                #out_string += "\t".join([chemical] + ['%d'%(rand_chem_count[score]) for score in chem_score_list[chemical]]) + '\n' 
#                out_string += chemical
#                for i in range(len(k_list)):
#                    paths_above_chem_score_count = rand_chem_count[chem_score_list[chemical][i]]
#                    out_string += "\t%0.6f\t%d" % (rand_k_score[k_list[i]], paths_above_chem_score_count)
#                out_string += "\n"
                out_string_med += chemical
                for k in k_list:
                    ties_count = sum([1 for curr_k in rand_k_score if rand_k_score[curr_k] >= rand_k_score[k]])
                    median_paths_score = np.median(list(rand_k_score.values())[:ties_count])
                    out_string_med += "\t%0.6f\t%d" % (median_paths_score, ties_count)
                out_string_med += '\n'

            #except (ValueError, IOError, IndexError):
            #    # if one of the chemicals failed, skip writing the output for this index so that it can be easily re-ran on the super computer
            #    print("ERROR: %s failed to get path stats. Not writing the file for this random_index." % (rand_paths))
            #    continue

            # write the output file
            out_file = rand_out_template % rand_index
            #print("writing %s" % (out_file))
            with open(out_file, 'w') as out:
                out.write(out_string_med)

    def proteinCountPval(self, chemicals_file, group_num_rectfs=False, inputs_dir=None, forced=False):
        if group_num_rectfs:
            print("--group-num-rectfs not yet setup for this option")
            sys.exit(1)
        print("Computing protein count pvalues for each chemical")
        self.chemicals = utils.readItemList(chemicals_file, col=1)
        t_utils.checkDir("%s" % self.out_dir)
        out_file = open("%s/chemical-protein-count-pvals.txt" % self.out_dir, 'w')
        # write the header line
        out_file.write("#chemical\tk\tnum_proteins\tpval\trand_num_proteins\n")

        for chemical in self.chemicals:
            print("\t%s parsing the protein counts" % (chemical),)

            #t_utils.checkDir("%s/%s" % (self.out_dir, chemical))
            #chem_out = "%s/%s/%s-protein-count-pvals.txt"%(self.out_dir, chemical, chemical)

            chem_ranked_edges = self.pathlinker_ranked_edges % chemical
            prot_count = self.getProteinCount(chem_ranked_edges, self.k_limit)
             
            print("\t%s parsing the rand protein counts" % (chemical),)
            rand_prot_counts = []
            # Get the random network's path probabilities
            for i in xrange(self.num_random):
                # TODO parallelize reading/getting the path probabilities from the random files 
                if group_num_rectfs:
                    rand_ranked_edges = self.random_pathlinker_ranked_edges %(num_rectfs, i+1)
                else:
                    rand_ranked_edges = self.random_pathlinker_ranked_edges %(chemical, i+1)
                rand_prot_count = self.getProteinCount(rand_ranked_edges, self.k_limit)
                rand_prot_counts.append(rand_prot_count) 

            #with open(chem_out, 'w') as out:
            #    out.write("#k\t

            # compute the p-val 
            prot_count_pval = t_utils.computePval(prot_count, rand_prot_counts, reverse=True)

            # write each chemicals output to the same file
            if self.num_random <= 100:
                formatted_pval = "%0.2f" % prot_count_pval
            elif self.num_random <= 1000:
                formatted_pval = "%0.3f" % prot_count_pval
            elif self.num_random <= 10000:
                formatted_pval = "%0.4f" % prot_count_pval

            out_file.write("%s\t%d\t%d\t%s\t%s\n" % (chemical, self.k_limit, prot_count, formatted_pval, ','.join([str(x) for x in sorted(rand_prot_counts, reverse=True)])))

        return

    def main(self, chemicals_file, group_num_rectfs=False, inputs_dir=None, forced=False):

        print("Computing pvalues for each chemical")
        self.chemicals = utils.readItemList(chemicals_file, col=1)

        # make sure the outptu directories are setup
        t_utils.checkDir("%s/viz"%self.out_dir)

        chem_pvals = self.out_dir + '/chemical-pvals.txt'
        # open these files and write the header line
        out_file = open(chem_pvals, 'w')
        out_file.write("#chemical\tX_k_pval (%d paths)\tmin_X_k\tmin_X_k_pval\tmin_k\tmin_k_pval\n"%self.k_limit)
        gpd_pval_out = open("%s/gpd_pvals.txt" % self.out_dir, 'w')
        gpd_pval_out.write("#chemical\tk\tprob\tperm pval\tgpd pval\tM\tnexc\tk\ta\tgood_fit_pval\texceedences\n")
        chem_rand_200_counts = open("%s/k200-counts.txt" % self.out_dir, 'w') 
        chem_rand_200_counts.write("#chemical\tk\tprob\tperm pval\trand_counts\n")

        if opts.group_num_rectfs:
            # get the # of rec and tfs perturbed by each chemical
            chem_num_rectfs = t_utils.build_chem_num_rectfs(self.chemicals, "%s/%s-rec-tfs.txt", inputs_dir)
            # get the chemicals grouped by the # of rec and tfs
            num_rectfs_chem = {} 
            for chem in chem_num_rectfs:
                if chem_num_rectfs[chem] not in num_rectfs_chem:
                    num_rectfs_chem[chem_num_rectfs[chem]] = set() 
                num_rectfs_chem[chem_num_rectfs[chem]].add(chem)
        else:
            # otherwise have each chemical map to itself
            num_rectfs_chem = {chem:[chem] for chem in self.chemicals}

        for num_rectfs in tqdm(sorted(num_rectfs_chem)):
            print(num_rectfs)
            # store the random results so they don't have to be read again
            loaded_random = False
            rand_k_probs = []
            rand_prob_counts = []
            #try:
            for chemical in num_rectfs_chem[num_rectfs]:

        ## Compute the significance of each chemical's network using path probability
        ##: Compute the fraction of random networks with greater path probability than the chemical's network
        #for chemical in tqdm(self.chemicals):

                t_utils.checkDir("%s/%s" % (self.out_dir, chemical))
                chem_out = "%s/%s/%s-pvals.txt"%(self.out_dir, chemical, chemical)
                # contains the # of paths at each prob cutoff, running total of paths, and the random median # of paths, 25 percentile, and 75 percentile
                prob_out = "%s/%s/%s-probs.txt"%(self.out_dir, chemical, chemical)
                # write each chemical's values to a file as well
                prob_counts_out = "%s/%s/%s-prob-counts.txt"%(self.out_dir, chemical, chemical)
                k_prob_out = "%s/%s/%s-k-prob.txt"%(self.out_dir, chemical, chemical)

                chem_k_prob = {}
                chem_prob_count = {}
                #if not os.path.isfile(chem_out) or forced:
                if os.path.isfile(prob_counts_out) and os.path.isfile(k_prob_out) and \
                        os.path.getsize(prob_counts_out) > 50 and os.path.getsize(k_prob_out) > 50 \
                        and not forced:
                    print("\t%s loading parsed chem and rand values" % (chemical),)
                    # load the already parsed chem and rand values
                    prob_counts_out_file = utils.readColumns(prob_counts_out, 1, 2, 3)
                    k_prob_out_file = utils.readColumns(k_prob_out, 1, 2, 3)
                    # the prob_counts file was written as cumulative counts, 
                    # so subtract the last amount from the current amount to get individual probability counts 
                    # which will be made into the cumulative counts later
                    last_count = 0
                    for prob, chem_count, rand_counts_test in prob_counts_out_file:
                        chem_prob_count[float(prob)] = int(chem_count) - last_count
                        last_count = int(chem_count)
                    for k, chem_prob, rand_probs in k_prob_out_file:
                        chem_k_prob[int(k)] = float(chem_prob)
                    #print(chem_k_prob)
                    if not loaded_random:
                        loaded_random = True 
                        last_counts = []
                        for i in xrange(self.num_random):
                            rand_k_probs.append({})
                            rand_prob_counts.append({})
                            last_counts.append(0)
                        for prob, chem_count, rand_counts_test in prob_counts_out_file:
                            rand_counts_test = rand_counts_test.split(',') 
                            for i in xrange(len(rand_counts_test)):
                                rand_prob_counts[i][float(prob)] = int(rand_counts_test[i]) - last_counts[i]
                                last_counts[i] = int(rand_counts_test[i])
                        for k, chem_prob, rand_probs in k_prob_out_file:
                            rand_probs = rand_probs.split(',') 
                            for i in xrange(len(rand_probs)):
                                rand_k_probs[i][int(k)] = float(rand_probs[i])
                else:
                    try:
                        print("\t%s parsing chem values" % (chemical),)
                        #print("Computing p-vals for %d %s"%(self.chemicals.index(chemical)+1, chemical))
                        #print("\tparsing files")
                        # first get the chemical's subnetwork path probabilities
                        chem_paths = self.pathlinker_paths % chemical
                        chem_k_prob, chem_prob_count = self.getPathStats(chem_paths, self.k_limit, self.group_by_prob)

                        if not loaded_random:
                            loaded_random = True 
                            rand_k_probs = []
                            rand_prob_counts = []
                            print("\t%s parsing rand values" % (chemical),)
                            # Get the random network's path probabilities
                            for i in xrange(self.num_random):
                                # TODO parallelize reading/getting the path probabilities from the random files 
                                if group_num_rectfs:
                                    rand_paths = self.random_pathlinker_paths %(num_rectfs, i+1)
                                else:
                                    rand_paths = self.random_pathlinker_paths %(chemical, i+1)
                                    if os.path.isfile(rand_paths) and os.path.getsize(rand_paths) < 50:
                                        print("Using %d in replace of problematic %d" % (i, i+1))
                                        # temporary fix to skip problematic samples
                                        rand_paths = self.random_pathlinker_paths %(chemical, i)
                                    if not os.path.isfile(rand_paths):
                                        #print("%s doesn't exist. Trying without the '-random'" % paths_file)
                                        # I accidentally wrote a bunch without the -random, so check to see if they still exist.
                                        new_rand_paths = rand_paths.replace('-random','')
                                        if os.path.isfile(new_rand_paths) and os.path.getsize(new_rand_paths) > 50:
                                            rand_paths = new_rand_paths
                                        else:
                                            print("Using %d in replace of problematic %d" % (i, i+1))
                                            # temporary fix to skip problematic samples
                                            rand_paths = self.random_pathlinker_paths %(chemical, i)
                                            if not os.path.isfile(rand_paths):
                                                rand_paths = rand_paths.replace('-random','')
                                rand_k_prob, rand_prob_count = self.getPathStats(rand_paths, self.k_limit, self.group_by_prob)
                                rand_k_probs.append(rand_k_prob)
                                rand_prob_counts.append(rand_prob_count)

                    except KeyError:
                        print("")
                        sys.stderr.write("\nFailed on chemical %s\n"%chemical)
                        raise

                #print("\tcomputing pvals")
                print("\t%s computing pvals" % (chemical),)
                if not self.group_by_prob: 
                    sum_pvals, k_pvals = self.computePathProbPval(chemical, chem_k_prob, rand_k_probs, self.k_limit, pval_sum=True)
                    # find the k value which has the minimum sum probability and minimum path probability
                    min_k_pval = (0,1)
                    min_sum_pval = (0,1)
                else:
                    sum_counts, rand_counts, prob_cum_count_pvals, prob_count_pvals = self.computePathProbPval(chemical, chem_prob_count, rand_prob_counts, self.k_limit, pval_sum=True, reverse_pval=self.group_by_prob)
                    # find the k value which has the minimum sum probability and minimum path probability
                    min_prob_pval = (0,1)
                    min_cum_pval = (0,1)
                print("\t%s writing output files" % (chemical),)
                # write the pvals at each k value for this chemical
                with open(chem_out, 'w') as out:
                    if not self.group_by_prob: 
                        for k in sorted(k_pvals):
                            out.write("%s\t%0.3f\t%0.3f\n" % (str(k)))
                            if k_pvals[k] < min_k_pval[1]:
                                min_k_pval = (k, k_pvals[k])
                            if sum_pvals[k] < min_sum_pval[1] and k > 25:
                                min_sum_pval = (k, sum_pvals[k])
                    else:
                        for prob in sorted(prob_count_pvals, reverse=True):
                            out.write("%0.6f\t%0.3f\t%0.3f\n" % (prob, prob_cum_count_pvals[prob], prob_count_pvals[prob]))
                            if prob_count_pvals[prob] < min_prob_pval[1]:
                                min_prob_pval = (prob, prob_count_pvals[prob])
                            # make sure the probability has at least 25 paths in it
                            #if prob_cum_count_pvals[prob] < min_cum_pval[1] and (prob > 25 or (self.group_by_prob and prob < 0.6)):
                            if prob_cum_count_pvals[prob] < min_cum_pval[1] and chem_k_prob[25] > prob:
                                min_cum_pval = (prob, prob_cum_count_pvals[prob])

                if self.group_by_prob:
                    self.write_probs(chem_prob_count, rand_prob_counts, prob_out)

                if self.group_by_prob:
                    # write the chemical, the very last (20,000 path) probability's pval, the minimum X_k pval, and the minimum pval
                    out_file.write("%s\t%0.6f\t%0.4f\t%0.6f\t%0.4f\t%0.4f\t%0.4f\n"%(chemical, sorted(prob_cum_count_pvals)[0], prob_cum_count_pvals[sorted(prob_cum_count_pvals)[0]], min_cum_pval[0], min_cum_pval[1], min_prob_pval[0], min_prob_pval[1]))
                else:
                    out_file.write("%s\t%0.3f\t%d\t%0.3f\t%d\t%0.3f\n"%(chemical, sum_pvals[self.k_limit], min_sum_pval[0], min_sum_pval[1], min_k_pval[0], min_k_pval[1]))

                # get the probability at a specific k, and then get all of the paths above that cutoff
                #k_prob = get_k_prob(chem_prob, k=200)
                #rand_prob_k = [rand_sums[i]+=rand_prob[i][k] for i in xrange(len(rand_prob))]
                chem_count = sum_counts[chem_k_prob[200]]
                #sum_counts, rand_counts,  
                #print("%s: prob: %s, count: %s, pval: %s" % (chemical, chem_k_prob[200], chem_count, prob_cum_count_pvals[chem_k_prob[200]]))
                rand_counts_prob = [rand_counts[chem_k_prob[200]][i] for i in xrange(len(rand_k_probs))]
                #print(rand_counts_prob)
                fig_file = "%s/viz/%s-k%dminsump-cum.png"%(self.out_dir, chemical, 200)
                #self.plotCumulativeHisto(chem_count, rand_counts_prob, fig_file, min_prob_pval[0], min_prob_pval[1])

                # write the chem_count and rand_count to use again later
                perm_pval = prob_cum_count_pvals[chem_k_prob[200]]
                chem_rand_200_counts.write("%s\t%d\t%0.4f\t%0.4f\t%s\n" % (chemical, chem_count, chem_k_prob[200], perm_pval, ','.join([str(x) for x in sorted(rand_counts_prob, reverse=True)])))

                # try the gdpPerm
                # Use 250 as the first Nmax to look at
                gpd_pval, M, k, a, nexc, good_fit_pval = gpdPerm.est(chem_count, rand_counts_prob, Nexcmax=250)
                # if M is >= 10, then the algorithm returns the p-value (permutation p-value) I already computed
                if M < 10:
                    print("\t%s computing the gpd p-value" % (chemical),)
                    gpd_pval_out.write("%s\t%d\t%0.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%0.4f\n" % (chemical, chem_count, chem_k_prob[200],
                            str(prob_cum_count_pvals[chem_k_prob[200]]), str(gpd_pval), str(M), str(nexc), str(k), str(a), good_fit_pval))
                # don't re-write the output files if they're already written
                if not os.path.isfile(prob_counts_out) or not os.path.isfile(k_prob_out) or \
                        os.path.getsize(prob_counts_out) < 50 or os.path.getsize(k_prob_out) < 50 \
                        or forced:
                    # write these values so they can be read again later if needed.
                    with open(prob_counts_out, 'w') as out:
                        out.write("#path score\tX_p count\tX_p rand_count\n")
                        for prob in sorted(prob_count_pvals, reverse=True):
                            # write the cumulative counts rather than the actual counts 
                            # because the random network cumulative counts have the same keys (probabilities) as the chemical. This makes it easier to write them all to a single file.
                            #out.write("%0.6f\t%d\t%s\n" % (prob, chem_prob_count[prob], ','.join([str(rand_probcount[prob]) for rand_probcount in rand_prob_counts])))
                            out.write("%0.6f\t%d\t%s\n" % (prob, sum_counts[prob], ','.join([str(count) for count in rand_counts[prob]])))
                    with open(k_prob_out, 'w') as out:
                        out.write("#path index\tpath score\trand path scores\n")
                        for k in chem_k_prob:
                            out.write("%d\t%0.6f\t%s\n" % (k, chem_k_prob[k], ','.join(["%0.6f"% rand_kprob[k] for rand_kprob in rand_k_probs])))

            #except IOError:
            #    # if the paths files are not completely written, skip for now
            #    sys.stderr.write("\nFailed on rec-tfs %s. Skipping for now\n" % (num_rectfs))
            #    continue 
            #         
            # make sure the python garbage collector cleans up unused memory after each loop
            gc.collect()

            #fig_file = "%s/viz/%s-%sminsump-cum.png"%(self.out_dir, chemical, str(min_sum_pval[0]))
            #self.plotCumulativeHisto(chem_sum, rand_sums, fig_file, min_sum_pval[0], min_sum_pval[1], sum_k=True)

            #print("\tplotting")
            # now do some plotting
            #chem_sum, rand_sums = self.getSums(chem_prob, rand_prob, min_sum_pval[0])
            ##rand_prob_k = [rand_sums[i]+=rand_prob[i][k] for i in xrange(len(rand_prob))]
            #fig_file = "%s/viz/%s-%sminsump-cum.png"%(self.out_dir, chemical, str(min_sum_pval[0]))
            #self.plotCumulativeHisto(chem_sum, rand_sums, fig_file, min_sum_pval[0], min_sum_pval[1], sum_k=True)
            # plot at an individual k value or probability value
            #values = [200, 20000] if not self.group_by_prob else [
            #for k in [200, 20000]:
            #    fig_file = "%s/viz/%s-%dcum.png"%(self.out_dir, chemical, k)
            #    rand_prob_k = [rand_prob[i][k] for i in xrange(len(rand_prob))]
            #    self.plotCumulativeHisto(chem_prob[k], rand_prob_k, fig_file, k, k_pvals[k])

        print("Done writing %s" % chem_pvals)
        out_file.close()
        gpd_pval_out.close() 
        chem_rand_200_counts.close() 

        # write the parsed random networks to a json file which can b e locded later
        #with open("%srectfs_rand_prob.json" % self.out_dir, 'w') as out_json:
        #    json.dump(rectfs_rand_prob, out_json, sort_keys=True, indent=2) 

        # Now print(some summary stats on the statistical significance)
        start_index = 2 if not self.group_by_prob else 3
        max_k_pvals = utils.readItemList(chem_pvals, start_index)
        sum_pvals = utils.readItemList(chem_pvals, start_index+2)
        k_pvals = utils.readItemList(chem_pvals, start_index+4)
        k200_pvals = utils.readItemList("%s/k200-counts.txt" % self.out_dir, 4)
        print
        print("At a p-value cutoff of %0.2f:" % self.pval_cutoff)
        print("%d chemical's paths are significant at k=%d for X_k pval"%(len([pval for pval in k200_pvals if float(pval) < self.pval_cutoff]), 200))
        print("%d chemical's paths are significant at k=%d for X_k pval"%(len([pval for pval in max_k_pvals if float(pval) < self.pval_cutoff]), self.k_limit))
        print("%d chemical's paths are significant at any k for X_k pvals"%(len([pval for pval in sum_pvals if float(pval) < self.pval_cutoff])))
        print("%d chemical's paths are significant at any k"%(len([pval for pval in k_pvals if float(pval) < self.pval_cutoff])))

    #def get_k_prob(self, chem_prob, k=200):
    #    for prob in sorted(chem_prob, reverse=True):

    def write_probs(self, chem_prob, rand_prob, out_file):
        """ Write the # of paths at each (probability) cutoff 
        Also write the median, and std deviation of the random networks
        """
        #print("writing %s" % out_file)
        # also write the X_k at each cutoff
        with open(out_file, 'w') as out:
            out.write("#prob_cutoff\t# of paths\trunning total\trandom_median\trandom_25th_percentile\trandom_75th_percentile\n")
            for c in sorted(chem_prob, reverse=True): 
                rand_sum_probs = [] 
                for i in xrange(len(rand_prob)):
                    # The path probabilities of the random networks might not match up with the chemical's network
                    # Compute the sum of the # of paths above the probability cutoff (k) (for example: # of paths above a path probability of 0.5)
                    rand_sum_probs.append(sum([rand_prob[i][prob] for prob in rand_prob[i] if prob >= c]))
                chem_sum = sum([chem_prob[prob] for prob in chem_prob if prob >= c])
                out.write("%0.10f\t%d\t%d\t%d\t%0.2f\t%0.2f\n" % (c, chem_prob[c], chem_sum, np.median(rand_sum_probs), np.percentile(rand_sum_probs, 25), np.percentile(rand_sum_probs, 75)))

    def getSums(self, chem_prob, rand_prob, k):
        # first get the k path probability
        chem_sum = 0 
        rand_sums = [0]*100
        for k in range(1,k+1):
            chem_sum += chem_prob[k] 
            for i in xrange(len(rand_prob)):
                rand_sums[i] += rand_prob[i][k]
        return chem_sum, rand_sums


if __name__ == "__main__":
    ## Parse command line args.
    usage = '%s [options]\n'%sys.argv[0]
    parser = OptionParser(usage=usage)
    parser.add_option('', '--chemicals', type='string', metavar='STR', default='results/chemicals/chemicals.txt',
                      help='Chemicals file. default=%default')
    parser.add_option('','--inputs-dir', type='string',
                      help='Dir containing the perturbed rec and tfs of each chemical. Required for --group-num-rectfs option.')   
    parser.add_option('','--forced',action='store_true', default=False,
                      help='Compute the statistical significance even if the files already exist')   
    #parser.add_option('','--group-num-rectfs',action='store_true', default=False,
    #                  help='Random network generation is grouped by the # of rec and tfs rather than random networks for each chemical.')   
    parser.add_option('', '--paths-dir', type='string', metavar='STR', default='results/chemicals/cyclinker/weight-no_change',
                      help='Dir of cyclinker results. Uses the paths file. default=%default')
    parser.add_option('', '--random-paths-dir', type='string', metavar='STR', default='results/random-sets/cyclinker',
                      help='Dir of random cyclinker results. Uses the paths file. default=%default')
#    parser.add_option('', '--out', type='string', metavar='STR', default='results/prec-rec.txt', \
#                      help='Output file')
    # these options are no long needed as pathlinker and cyclinker output the same file format
#    parser.add_option('-P', '--pathlinker', action='store_true', default=False,
#                      help='Paths are in pathlinker format')
#    parser.add_option('-C', '--cyclinker', action='store_true', default=False,
#                      help='Paths are in cyclinker format')
    parser.add_option('', '--k-limit', type='int', default=200,
                      help='limit of k paths for computations. default=%default')
    parser.add_option('-k', '--k-to-test', type='int', action='append', 
                      help='Value of k to test for significance. Multiple k values can be given.')
    parser.add_option('', '--group-by-prob', action='store_true', default=False,
                      help='group the p-values by probability values instead of k values ')
    parser.add_option('-n', '--num-random', type='int', default=('1','100'),nargs=2, 
                      help='The number of random sets sets to use. Default= 1 100')
    parser.add_option('-o', '--out-dir', action='store',
                      help='Output dir to place the files')
    parser.add_option('', '--protein-count-pval', action='store_true', default=False,
                      help='Compute a p-val for the number of proteins or nodes in the top --k-limit paths')
    parser.add_option('', '--write-counts', action='store_true', default=False,
                      help='Write the path score at each given k for each chemical')
    parser.add_option('', '--write-rand-counts', action='store_true', default=False,
                      help='Write the number of paths in each random network at each chemicals score')
    parser.add_option('', '--compute-score-pval', action='store_true', default=False,
                      help='Parse the chemical and random counts and compute a permutation and gpd pval for each specified k')
    parser.add_option('', '--use-counts', action='store_true', default=False,
                      help='Parse the chemical and random counts and compute a permutation and gpd pval for each specified k')

    (opts, args) = parser.parse_args()

    stat_sig = StatSig(opts.paths_dir, opts.out_dir, opts.random_paths_dir, opts.k_limit, opts.num_random, opts.group_by_prob)
    if opts.protein_count_pval:
        stat_sig.proteinCountPval(opts.chemicals, opts.group_num_rectfs, opts.inputs_dir, opts.forced)
    if opts.write_counts:
        stat_sig.write_counts(opts.chemicals, opts.k_to_test, opts.forced)
    if opts.write_rand_counts:
        stat_sig.write_rand_counts(chemicals_file=opts.chemicals, forced=opts.forced)
    if opts.compute_score_pval:
        stat_sig.compute_score_pval(opts.k_to_test, opts.use_counts, opts.forced)
    # just run the main if none of the other options were specified
    if not opts.protein_count_pval and not opts.write_counts and not opts.write_rand_counts and not opts.compute_score_pval:

#        if (not opts.cyclinker and not opts.pathlinker) or (opts.cyclinker and opts.pathlinker):
#            print("Error: must specify either --cyclinker or --pathlinker option")
#            parser.print_help()
#            sys.exit(1)

        if opts.group_num_rectfs and not opts.inputs_dir:
            print("Error: --inputs-dir is required for the --group-num-rectfs option")
            sys.exit(1)
        
        stat_sig.main(opts.chemicals, opts.group_num_rectfs, opts.inputs_dir, opts.forced)

