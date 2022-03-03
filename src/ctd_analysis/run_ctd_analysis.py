
print("Importing libraries")
from collections import defaultdict
from optparse import OptionParser,OptionGroup
import os
import sys
import pandas as pd
from tqdm import tqdm
from src import toxcast_utils as t_utils
from scipy.stats import hypergeom
import gzip



def parse_args():
    usage = 'run_ctd_analysis.py [options]\n'
    # 
    # run the DAVID analysis on each of them
    # and write each to a file, as well as one combined file(?)
    description = "Compute the statistical significance of the overlap of the the proteins in our networks and the proteins in involved in phosphorylation reactions in CTD as a result of chemical exposure"
    parser = OptionParser(usage=usage)

    ## Common options
    parser.add_option('','--chemicals',type='string',
                      help='File containing chemicals to analyze. Required')
    parser.add_option('','--ctd-file',type='string',
                      help='File containing the CTD reactions. Required')
    parser.add_option('','--interactome',type='string',
                      help='Interactome file from which to get the universe of proteins. Required')
    parser.add_option('','--mapping-file',type='string',
                      help='Uniprot to gene name mapping file')
    parser.add_option('-S','--single-run', type='string', action='append',
            help='Run only a single chemical. Can specify multiple chemicals by listing this option multiple times.')
    parser.add_option('-o','--out-dir', action='store',
                      help='output to place the output files in. Required')
    parser.add_option('','--paths-dir',type='string',metavar='STR',
            help='Dir of Cyclinker results. Uses the paths file. Example: outputs/<version>/weighted/edgelinker/. Required')
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
    if (not opts.chemicals and not opts.single_run) or not opts.paths_dir \
       or not opts.out_dir or not opts.ctd_file: 
        print("Must specify --chemicals, --paths-dir, --out-dir, --ctd-file. Exiting")
        parser.print_help()
        sys.exit(1)

    return opts

# will be loaded in main
toxcast_data = None
uniprot_to_gene = None

def main(chemicals, paths_dir, out_dir, 
         pval_cutoff=0.05, corr_type="BF", **kwargs):

    global toxcast_data, uniprot_to_gene
    toxcast_data = t_utils.loadToxcastData()
    chemIDtoCAS, chemCAStoID = get_chemical_map(toxcast_data)
    uniprot_to_gene_df = pd.read_csv(kwargs['mapping_file'], sep='\t', header=None)
    uniprot_to_gene = dict(uniprot_to_gene_df.values)
    t_utils.checkDir(out_dir)

    chem_prots, reports_dir = load_prots(chemicals, paths_dir, out_dir, **kwargs)

    ctd_genes, ctd_chem_itxs = load_ctd_data(kwargs['ctd_file'], chemCAStoID)

    # To get the background set of genes for the hypergeometric test, 
    # get the proteins that are both in CTD and in the interactome
    print("reading %s" % (kwargs['interactome']))
    df = pd.read_csv(kwargs['interactome'], sep='\t',comment='#', header=None)
    ppi_prots = set(df[0]) | set(df[1])
    # map both to the gene name space
    ppi_genes = set(uniprot_to_gene[p] for p in ppi_prots)
    print("\t%d interactome_genes" % (len(ppi_genes)))

    background_genes = ppi_genes & ctd_genes
    print("%d genes both in the interactome and in CTD" % (len(background_genes)))
    print("limiting CTD phosphorylation interactions to those in the interactome")
    ctd_chem_itxs = {c: p & ppi_genes for c, p in ctd_chem_itxs.items()}
    # also add the other chemicals
    for c in chemicals:
        if c not in ctd_chem_itxs:
            ctd_chem_itxs[c] = set()

    chem_pval = {}
    chem_net_prots_with_ctd = {} 
    # TODO try both making random subsets and the hypergeometric test
    pop_size = len(background_genes)
    #num_success_states_in_pop = len(set(p for c,p in ctd_itxs.items()))
    for chem, prots in chem_prots.items():
        genes = set(uniprot_to_gene[p] for p in prots)
        #if len(genes) != len(prots):
        #    print("Warning: %s: num genes != num prots! (%s, %s)" % (chem, len(genes), len(prots)))
        num_genes_with_ctd = len(genes & ctd_chem_itxs[chem])
        chem_net_prots_with_ctd[chem] = num_genes_with_ctd
        # number of draws is the # genes in the network
        num_draws = len(genes)
        # number of successes is the # genes in the network with a CTD interaction
        num_successes = num_genes_with_ctd
        # number of success stats in the population is the number of phosphorylation interactions of this chemical
        num_success_states_in_pop = len(ctd_chem_itxs[chem])
        M, n, N, k = pop_size, num_success_states_in_pop, num_draws, num_successes
        # Use k-1 since the survival function (sf) gives 1-cdf. The cdf at k gives the probability of drawing k or fewer. The sf at k is the probability of drawing k+1 or more
        # https://blog.alexlenail.me/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
        # https://github.com/scipy/scipy/issues/7837
        pval = hypergeom.sf(k-1, M, n, N)
        chem_pval[chem] = pval

    # now write to a file
    out_file = "%s/CTD-stat-sig.tsv" % (out_dir)
    print("writing %s" % (out_file))
    with open(out_file, 'w') as out:
        header_line = '\t'.join([
            "ChemID", "ChemName", "# net prots", "# CTD phospho prots",
            "# overlap", "pval", "BF corr-pval"])
        out.write(header_line + '\n')
        for chem, prots in chem_prots.items():
            name = toxcast_data.chemIDtoName[chem]
            out.write('\t'.join(str(x) for x in [
                chem, name, len(prots), len(ctd_chem_itxs[chem]),
                chem_net_prots_with_ctd[chem], chem_pval[chem],
                chem_pval[chem] * len(chemicals)]) + '\n')


def get_chemical_map(toxcast_data):
    """
    Map the chemical ID used in ToxCast to the CAS RN 
    """
    print("reading %s" % (toxcast_data.chemical_summary_file))
    df = pd.read_csv(toxcast_data.chemical_summary_file, header=0)
    cols = df.columns
    chemIDtoCAS  = dict(zip(df[cols[3]], df[cols[2]]))
    chemCAStoID  = dict(zip(df[cols[2]], df[cols[3]]))
    return chemIDtoCAS, chemCAStoID


def load_ctd_data(ctd_file, chemCAStoID):
    # read the ctd interactions
    ctd_itxs = defaultdict(set)
    ctd_genes = set()
    print("reading %s" % (ctd_file))
    with gzip.open(ctd_file, 'r') as f:
        for line in f:
            line = line.decode() 
            if line[0] == '#':
                continue
            line = line.rstrip().split('\t')
            chem_cas_id = line[2]
            # map these to the other chemical namespace we use for toxcast if available
            chem_id = chemCAStoID.get(chem_cas_id, chem_cas_id)
            # just use the gene name for now(?)
            gene = line[3]
            gene_form = line[5]
            taxon = line[7]
            interaction = line[8]
            if taxon != "9606" or gene_form != "protein": 
                continue
            ctd_genes.add(gene)
            if 'phosphorylation' in interaction:
                ctd_itxs[chem_id].add(gene) 
    print("\t%d unique chemical-protein phosphorylation interactions" % (len([itx for itx in ctd_itxs.items()])))
    return ctd_genes, ctd_itxs


def load_prots(chemicals, paths_dir, out_dir, k_limit=200, **kwargs):
    if kwargs['run_on_hits']:
        chem_prots = toxcast_data.chemical_protein_hit
        reports_dir = "%s/chemical-hits-reports/" % (out_dir)
    else:
        # one chemical and protein ID pair on each line
        prots_file = "%s/chem-prots.txt" % (out_dir)
        if not kwargs['forced'] and os.path.isfile(prots_file):
            print("reading %s. Use --forced to overwrite" % (prots_file))
            s = pd.read_csv(prots_file, sep='\t', index_col=0, header=None, squeeze=True)
            # now convert it back to a dictionary
            chem_prots = {chem: prots.to_list() for chem, prots in s.groupby(s.index)}
        else:
            # load the proteins in each chemical's network
            edgelinker_output = paths_dir+'/%s-paths.txt'
            print("Reading paths for each chemical network from: %s" % (edgelinker_output))
            chem_prots = {}
            for chemical in chemicals:
                proteins = t_utils.getProteins(paths=edgelinker_output % chemical, max_k=k_limit)
                chem_prots[chemical] = list(proteins)
            s = pd.Series(chem_prots).explode()
            print("writing %s" % (prots_file))
            s.to_csv(prots_file, sep='\t', header=False)

        reports_dir = "%s/chemical-reports/" % (out_dir)
    t_utils.checkDir(os.path.dirname(reports_dir))
    return chem_prots, reports_dir


if __name__ == "__main__":
    opts = parse_args()
    kwargs = vars(opts)
    if opts.single_run:
        chemicals = opts.single_run
    else:
        chemicals = pd.read_csv(opts.chemicals, sep='\t', header=None, comment='#')[0].tolist()
    del kwargs['chemicals']
    main(chemicals, **kwargs)
