
print("Importing libraries")
from optparse import OptionParser,OptionGroup
import os
import sys
import pandas as pd
from tqdm import tqdm
from src import toxcast_utils as t_utils
from src.goterm_analysis import david_client


def parse_args():
    usage = 'run_david_analysis.py [options]\n'
    # 
    # run the DAVID analysis on each of them
    # and write each to a file, as well as one combined file(?)
    description = "load the proteins in each chemical's network, run the DAVID analysis (GO_BP_DIRECT) on each of them, then analyze the overlap of the significant terms"
    parser = OptionParser(usage=usage)

    ## Common options
    parser.add_option('','--chemicals',type='string',
                      help='File containing chemicals to analyze. Required')
    parser.add_option('-S','--single-run', type='string', action='append',
            help='Run only a single chemical. Can specify multiple chemicals by listing this option multiple times.')
    parser.add_option('-o','--out-dir', action='store',
                      help='output to place the output files in. Required')
    parser.add_option('','--paths-dir',type='string',metavar='STR',
            help='Dir of Cyclinker results. Uses the paths file. Example: outputs/<version>/weighted/cyclinker/. Required')
    parser.add_option('-k','--k-limit',action='store', type='int', default=200,
                      help='k limit for the # paths to read. Default: 200')
    parser.add_option('','--run-on-hits',action='store_true', default=False,
                      help='Instead of the prots in the network, run a GO term analysis on the chemical hits.')
    #parser.add_option('', '--pdf', action='store_true',
    #                  help='Also store a pdf of the figures')
    parser.add_option('','--forced',action='store_true',default=False,
                      help='Re-load the proteins and re-run the DAVID analysis if its already been run')

    parser.add_option('','--pval-cutoff',type='float', default=0.05,
                      help="Cutoff on the significance of the enrichment of that term.")
    parser.add_option('','--correction-type',type='string', 
            help="p-value correction type. Options are: BF (Bonferroni), BH (Benjamini-Hochberg). Default: None")

    # parse the command line arguments
    (opts, args) = parser.parse_args()

    #if not opts.chemicals and not opts.sig_cutoff:
    if (not opts.chemicals and not opts.single_run) or not opts.paths_dir or not opts.out_dir: 
        print("Must specify --chemicals, --paths-dir, and --out-dir. Exiting")
        parser.print_help()
        sys.exit(1)

    return opts


def main(chemicals, paths_dir, out_dir, k_limit=200, 
         forced=False, pval_cutoff=0.05, corr_type="BF", **kwargs):

    t_utils.checkDir(out_dir)

    # one chemical and protein ID pair on each line
    prots_file = "%s/chem-prots.txt" % (out_dir)

    if kwargs['run_on_hits']:
        toxcast_data = t_utils.loadToxcastData()
        chem_prots = toxcast_data.chemical_protein_hit
        # limit the prots to those that are hit
        chem_prots = {c: [p for p, hit in hit_dict.items() if hit == 1] 
                      for c, hit_dict in chem_prots.items()}
        reports_dir = "%s/chemical-hits-reports/" % (out_dir)
    else:
        if not forced and os.path.isfile(prots_file):
            print("reading %s. Use --forced to overwrite" % (prots_file))
            s = pd.read_csv(prots_file, sep='\t', index_col=0, header=None, squeeze=True)
            # now convert it back to a dictionary
            chem_prots = {chem: prots.to_list() for chem, prots in s.groupby(s.index)}
        else:
            # load the proteins in each chemical's network
            cyclinker_output = paths_dir+'/%s-paths.txt'
            print("Reading paths for each chemical network from: %s" % (cyclinker_output))
            chem_prots = {}
            for chemical in chemicals:
                proteins = t_utils.getProteins(paths=cyclinker_output % chemical, max_k=k_limit)
                chem_prots[chemical] = list(proteins)
            s = pd.Series(chem_prots).explode()
            print("writing %s" % (prots_file))
            s.to_csv(prots_file, sep='\t', header=False)

        reports_dir = "%s/chemical-reports/" % (out_dir)
    t_utils.checkDir(os.path.dirname(reports_dir))
    # run the DAVID analysis on each of them
    client = None
    # reset the client each time. Maybe this isn't needed?
    for chem in tqdm(chemicals):
        chart_file = "%s/%s.txt" % (reports_dir, chem)
        if not forced and os.path.isfile(chart_file):
            print("%s already exists. Use --forced to overwrite" % (chart_file))
            continue
        if client is None:
            print("Setting up david client")
            client = david_client.DAVIDClient()
            client.set_category('GOTERM_BP_DIRECT')
        print(chem)
        prots = chem_prots[chem]
        # pass the list of proteins 
        client.setup_inputs(','.join(prots), idType='UNIPROT_ACCESSION', listName=chem)
        # make sure we're using the right list
        #print(client.client.service.getCurrentList())
        # build the functional annotation chart
        client.build_functional_ann_chart()

        # and write each to a file
        #print("writing %s" % (chart_file))
        client.write_functional_ann_chart(chart_file)

    pval_col = "Pvalue"
    if corr_type == "BF":
        pval_col = "Bonferroni"
    elif corr_type == "BH":
        pval_col = "Benjamini"
    # now read each of them and write a combined file
    dfs = []
    for chem in chemicals:
        chart_file = "%s/%s.txt" % (reports_dir, chem)
        if not os.path.isfile(chart_file):
            print("%s doesn't exist. Skipping" % (chart_file))
            continue
        df = pd.read_csv(chart_file, sep='\t')
        # apply the p-value cutoff
        df = df[df[pval_col] < pval_cutoff]
        df = df[['Term', pval_col]]
        # split the name and id
        df['GOID'] = df['Term'].apply(lambda x: x.split('~')[0])
        df['Term'] = df['Term'].apply(lambda x: x.split('~')[1])
        df['Chemical'] = chem
        print(len(df))
        dfs.append(df)

    df_all = pd.concat(dfs)
    print(df_all.head())
    all_terms_file = "%s/%schemical%s-sig-terms-%s-c%s.tsv" % (
        out_dir, len(chemicals), "-hits" if kwargs['run_on_hits'] else "s",
        pval_col.lower(), str(pval_cutoff).replace('.','_'))
    df_all.to_csv(all_terms_file, sep='\t', index=None,
                  columns=['Chemical', 'Term', 'GOID', pval_col])

    # now compare the overlap of the enriched terms!
    #df_all.groupby('Term').value_counts()
    counts = df_all[['Term','GOID']].value_counts()
    print(counts)
    counts_file = "%s/%schemicals-sig-terms-%s-c%s-counts.tsv" % (
        out_dir, len(chemicals), pval_col.lower(), str(pval_cutoff).replace('.','_'))
    print("writing to %s" % (counts_file))
    counts.to_csv(counts_file, header=False, sep='\t')

    # now run REVIGO on this list of terms?


if __name__ == "__main__":
    opts = parse_args()
    kwargs = vars(opts)
    if opts.single_run:
        chemicals = opts.single_run
    else:
        chemicals = pd.read_csv(opts.chemicals, sep='\t', header=None)[0].tolist()
        chemicals = [c for c in chemicals if '#' not in c]
    del kwargs['chemicals']
    main(chemicals, **kwargs)
