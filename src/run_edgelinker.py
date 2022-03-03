#!/usr/bin/env python
""" Script for running Peter's java implementation of edgelinker.
These functions were pulled out into their own script so that instances of edgelinker could be run in parallel.
"""

from optparse import OptionParser
import os
import sys
import subprocess
from src import toxcast_utils as t_utils


# out_pref is the path/to/output/name that will be used to name the *-paths.txt and *-ranked-edges.txt files
def runEdgeLinker(interactome, rec_tfs_file, out_pref, max_k=1000,
        multi_run=False, edge_penalty=None, rec_tfs_penalty=False):
    """
    *max_k*: the maximum amount of paths to let edgelinker write. I arbitrarily chose 60,000 as the default to give a little more space than the usual 20,000 default
    """
    #print("running edgelinker on %s " % edgelinker_paths_file)
    #interactome = os.path.abspath(interactome)
    #rec_tfs_file = os.path.abspath(rec_tfs_file)
    #out_pref = os.path.abspath(out_pref)
    #print(("Runnin EdgeLinker on %s"%out_pref))

    # now run edgelinker
    class_path = "src/EdgeLinker/bin:src/EdgeLinker/bin/commons-cli-1.3.1.jar"
    cmd = "java -cp %s EdgeLinker --network %s --start-ends %s --out-pref %s --max-k %d" % (class_path, interactome, rec_tfs_file, out_pref, max_k)
    if multi_run is True:
        cmd += " --multi-run "
    if edge_penalty is not None:
        cmd += " --edge-penalty %s " % (str(edge_penalty))
    if rec_tfs_penalty is True:
        cmd += " --start-ends-penalty "
    t_utils.runCommand(cmd, error_message="ERROR: EdgeLinker failed on %s!\n"%out_pref)

    # if there is nothing in the output file besides the header line, print(the failure as well)
    if not multi_run and os.path.getsize("%s-paths.txt"%out_pref) < 50:
        sys.stderr.write(cmd + '\n')
        sys.stderr.write("EdgeLinker failed on %s!\n"%out_pref)


# this function will utilize the edgelinker functionality of reading the interactome once
# (which takes about 2 seconds each time) and then running edgelinker on each of the chemicals
def createGroupEdgeLinker(chemicals, inputs_dir, out_pref):
    """
    *out_pref*: will create infiles.txt and outfiles.txt using the output prefix
    which contain all of the rec-tfs input files and the output prefix for each chemical's edgelinker results
    """
    # now run edgelinker on each of the chemicals using the permuted network
    rec_tfs_file_template = "%s/rec-tfs/%%s-rec-tfs.txt" % (inputs_dir) 
    in_files = []
    out_files = []
    # add all input files and output prefixes into a list
    for chemical in sorted(chemicals):
        rec_tfs_file = rec_tfs_file_template % (chemical)
        in_files.append(os.path.abspath(rec_tfs_file))
        out_files.append(os.path.abspath("%s%s" % (out_pref, chemical)))

    # write the in and out files to the networks dir
    edgelinker_in_files = '%sinfiles.txt' % (out_pref)
    with open(edgelinker_in_files, 'w') as out:
        out.write('\n'.join(in_files))
    edgelinker_out_files = '%soutfiles.txt' % (out_pref)
    with open(edgelinker_out_files, 'w') as out:
        out.write('\n'.join(out_files))
    return edgelinker_in_files, edgelinker_out_files


if __name__ == "__main__":
    ## Parse command line args.
    usage = '%s [options]\n'%sys.argv[0]
    parser = OptionParser(usage=usage)
    parser.add_option('', '--network', 
            help="Interactome to use when running edgelinker")
    parser.add_option('', '--rec-tfs-file', 
            help="File containing the receptors and tfs to use as sources and targets when running edgelinker")
    parser.add_option('', '--out-pref', 
            help="path/to/prefix for the output files.")
    parser.add_option('','--max-k',type='int',default='1000',
            help="maximum number of paths to keep when running edgelinker. Useful to not generate too much data (~2TB if all paths of each random network are kept)\n\tDefault = 1000")
    parser.add_option('', '--multi-run', action='store_true',
            help="Option to create two files containing the rec-tfs file and output-prefixes for all of the specified chemicals")
    parser.add_option('', '--edge-penalty', type='float', default='1',
            help="Add the natural log of the specified penalty to the cost of each edge. This will effectively increase the cost of each path by the length * edge penalty")
    parser.add_option('', '--rec-tfs-penalty', action='store_true',
            help="Set the cost specified in the third column of the rec-tfs-file to be the cost of the super-source->start and end->super-target edges. Otherwise will be 0")

    (opts, args) = parser.parse_args()

    if opts.network is None or opts.rec_tfs_file is None or opts.out_pref is None:
        print("\t--network, --rec-tfs-file and --out-pref are all required")
        sys.exit()

    runEdgeLinker(opts.network, opts.rec_tfs_file, opts.out_pref, opts.max_k, 
            multi_run=opts.multi_run, edge_penalty=opts.edge_penalty, rec_tfs_penalty=opts.rec_tfs_penalty)

    # TODO add actual command line arguments
#    if len(sys.argv) < 5:
#        print("Usage: python runEdgeLinker.py <interactome> <rec_tfs_file> <out_pref> <max_k>")
#        sys.exit()
#    if len(sys.argv) > 5 and sys.argv[5] == '--multi-run':
#        runEdgeLinker(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), multi_run=True)
#    else:
#        runEdgeLinker(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))

