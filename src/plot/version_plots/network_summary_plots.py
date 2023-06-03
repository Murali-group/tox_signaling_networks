# Quick script to plot figures about the proteins in the response networks

print("Importing modules")

import sys
import os
import shutil
from optparse import OptionParser
from src.utils.file_utils import readItemList
from src import toxcast_utils as t_utils
from src import toxcast_settings as t_settings
from src import summary_stats
# Plotting imports:
# To save files remotely, use Agg.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# try the ggplot style
#plt.style.use('ggplot')
import numpy as np
import seaborn as sns
sns.set_style("whitegrid")


def parseArguments():
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    ## Common options
    parser.add_option('','--version',type='string',action="append",
            help="Version of the PPI to run. Multiple can be run. Options are: %s." % (', '.join(t_settings.ALLOWEDVERSIONS)))
    parser.add_option('-c', '--compare-versions', action='store_true',
            help='plots will be written to viz/version_plots/summary-stats/summary-network-stats-version.png instead of the ' + 
            'default outputs/version/weighted/plots/summary-stats/summary-network-stats-version.png. If not specified, default is used.')
    parser.add_option('-o', '--out-file', type='string', metavar='STR',
            help='path/to/output_file.png. If not specified, default is used.')
    parser.add_option('-k', '--k-to-test', type='int',
            help='K to use to get the signaling networks')
    parser.add_option('', '--pdf', action='store_true',
            help='Also store a pdf of the figures')
    parser.add_option('','--poster',action="store_true",
            help="Create the version of this figure for the poster")
    parser.add_option('','--forced',action="store_true",
            help="Overwrite summary stats and figure if they already exist")

    # parse the command line arguments
    (opts, args) = parser.parse_args()
    return opts


opts = parseArguments()

#if __name__ != "__main__":
#    os.chdir("/data/jeff-law/projects/2016-02-toxcast/")

for version in opts.version:
    print("")
    print("-"*30)
    t_settings.set_version(version)
    chemicals = sorted(readItemList("%s/chemicals.txt" % (t_settings.INPUTSPREFIX)))
    #sig_chemicals = utils.readItemList("inputs/%s/sig-chemicals.txt" % (version), 1)
    #unsig_chemicals = set(chemicals).difference(set(sig_chemicals))

    #summary_file = "outputs/%s/weighted/stats/network-summaries.csv" % (version)
    summary_file = "network_summaries_k%s.csv"% (opts.k_to_test)
    df = summary_stats.get_summary_stats(
            version=version, 
            summary_file=summary_file,
            k_to_test=opts.k_to_test,
            forced=opts.forced)
    #df = t_utils.get_summary_stats(version=version, forced=True)
    print(df.drop(columns=["net_hits", "net_nonhits", "inter_net_hits", "inter_net_nonhits", "inter_hits", "inter_nonhits", "pvals", "qvals"]).describe())
    print(df[["net_hits", "net_nonhits", "inter_net_hits", "inter_net_nonhits", "inter_hits", "inter_nonhits", "pvals", "qvals"]].describe())

    # loop through the chemicals, significant chemicals and unsignificant chemicals
    for chemicals, postfix in [(chemicals, '')]:
        if opts.out_file is None:
            out_file_name = "summary-network-stats-%s-k%s.png" % (version, opts.k_to_test)
            out_dir = "%s/plots/summary-stats/" % (t_settings.RESULTSPREFIX)
            t_utils.checkDir(out_dir)
            out_file = "%s/%s" % (out_dir, out_file_name)
            # if specified, copy the file to the compare versions dir
            if opts.compare_versions:
                out_dir_compare_versions = "viz/version_plots/summary-stats"
                t_utils.checkDir(out_dir_compare_versions)
                out_file_compare_versions = "%s/%s" % (out_dir_compare_versions, out_file_name)
        else:
            out_file = opts.out_file

        print("Writing %s" % (out_file))
        #fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        # set the ratio for each of the sub figures
        if opts.poster is True:
            fig = plt.figure(figsize=(12, 5))
        else:
            fig = plt.figure(figsize=(12.5, 5))
        color_palette = sns.color_palette()
        gs = matplotlib.gridspec.GridSpec(1, 4, width_ratios=[4, 3, 1, 2])
        #fig.suptitle("Chemical Response Network Overview", fontsize=16)
        fig.text(0.03, 0.97, "A", horizontalalignment='left', verticalalignment='center', fontsize=16, weight="bold")
        fig.text(0.38, 0.97, "B", horizontalalignment='left', verticalalignment='center', fontsize=16, weight="bold")
        fig.text(0.66, 0.97, "C", horizontalalignment='left', verticalalignment='center', fontsize=16, weight="bold")
        fig.text(0.79, 0.97, "D", horizontalalignment='left', verticalalignment='center', fontsize=16, weight="bold")
        #plt.title("Chemical Hits and Response Network Overview")
        #plt.ylabel("Counts")
        df = df.rename(columns={
            'hit_rec': 'Responsive\nReceptors', 'hit_tfs': 'Responsive\nTFs',
            'prots': 'Network\nProteins', 'hits': 'Responsive\nProteins', 'nonhits': 'Non-Responsive\nProteins',
            'num_edges': "Network\nEdges", 'num_paths': "Network\nPaths", "avg_path_lengths": "Average\nPath Length",
            })
        df['Responsive\nRatios'] = df['inter_net_hits'] / df['inter_hits']
        df['Non-Responsive\nRatios'] = df['inter_net_nonhits'] / df['inter_nonhits']
        df = df.replace(np.inf,0).replace(np.nan,0)
        #plt.boxplot([num_prots, num_net_rec, num_net_tfs], labels=['Proteins', 'Receptors', 'TFs'])
        #ax1.boxplot([num_prots, num_hits, num_nonhits], labels=['Proteins', 'Hits', 'Non-Hits'])
        ax1 = plt.subplot(gs[0])
        sns.boxplot(data=df[["Responsive\nReceptors", "Responsive\nTFs", "Responsive\nProteins", "Non-Responsive\nProteins"]], ax=ax1)
        ax1.set_ylabel("Frequency")
        boxes = ax1.artists
        boxes[2].set_facecolor(color_palette[2])
        boxes[3].set_facecolor(color_palette[3])
        #y0, y1 = ax1.get_ylim()
        #ax1.set_ylim(0, y1)
        ax2 = plt.subplot(gs[1])
        sns.boxplot(data=df[["Network\nProteins", "Network\nPaths", "Network\nEdges"]], ax=ax2, palette="Set2")
        #ax2.boxplot([num_rec, num_tfs], labels=['Hit Receptors', 'Hit TFs'])
        ax2.set_ylabel("Frequency")
        #boxes = ax2.artists
        #boxes[2].set_facecolor(color_palette[4])
        ax3 = plt.subplot(gs[2])
        sns.boxplot(data=df[["Average\nPath Length"]], ax=ax3)
        ax3.set_ylabel("Frequency")
        ax4 = plt.subplot(gs[3])
        sns.boxplot(data=df[['Responsive\nRatios', 'Non-Responsive\nRatios']], ax=ax4)
        ax4.set_ylabel("Ratios")
        ax4.set_ylim(-0.02, 0.6)
        boxes = ax4.artists
        boxes[0].set_facecolor(color_palette[2])
        boxes[1].set_facecolor(color_palette[3])
        #plt.xlabel("In the Response Networks")
        plt.tight_layout()
        if opts.forced or not os.path.isfile(out_file):
            plt.savefig(out_file)
        else:
            print("figure already exists. use --forced to overwrite it")
        if opts.pdf and (opts.forced or not os.path.isfile(out_file.replace('.png', '.pdf'))):
            plt.savefig(out_file.replace('.png', '.pdf'))
            plt.savefig(out_file.replace('.png', '.svg'))
        plt.close()

        # print the stats for the case study chemicals
        measures = ["Responsive\nReceptors", "Responsive\nTFs", "Network\nProteins", 
                "Responsive\nProteins", "Non-Responsive\nProteins", "Responsive\nRatios", "Non-Responsive\nRatios",
                "Average\nPath Length", "Network\nPaths", "Network\nEdges", "qvals",
                ]
        # Lovastatin, T3, BPA, Cyclopamine
        #case_studies = ['23216', '20182', '20784']
        case_studies = ["C75330755", "C6893023", "C80057", "C4449518"]
        df = df[measures].loc[case_studies]
        #out_dir = "outputs/%s/weighted/stats/summary" % (version) 
        out_dir_csv = "%s/stats/summary" % (t_settings.RESULTSPREFIX) 
        t_utils.checkDir(out_dir_csv)
        df.to_csv("%s/case_study_summaries.csv" % (out_dir_csv))
    #    # --------------------------------------------------
        if opts.compare_versions and opts.out_file is None:
            out_file2 = out_file_compare_versions
            print("Copying '%s' to '%s'" % (out_file, out_file2))
            shutil.copy(out_file, out_file2)
            if opts.pdf:
                print("Copying '%s' to '%s'" % (out_file.replace('.png','.pdf'), out_file2.replace('.png','.pdf')))
                shutil.copy(out_file.replace('.png','.pdf'), out_file2.replace('.png','.pdf'))
                shutil.copy(out_file.replace('.png','.svg'), out_file2.replace('.png','.svg'))
