#! /usr/bin/env python
""" This script parses the chemical, assay, and hits matrix to get the hit/nonhit proteins as well as rec/tfs for each chemical. 
Is imported and used by many other toxcast scripts (like the master-script.py). 
Also writes hits files as well as the local and global random rec-tfs inputs files
"""

from collections import defaultdict
#import src.toxcast_utils as t_utils
from src import toxcast_utils as t_utils
import sys
from optparse import OptionParser
import os
import json
import pandas as pd
import numpy as np
from tqdm import tqdm
#import csbdb


class ToxCastData:
    #def __init__(self, interactome="pathlinker-data/2015pathlinker-weighted.txt", include_nuclear_receptors=False):
    def __init__(self, include_nuclear_receptors=False, forced=False, verbose=False):
        # boolean value to either include or exclude nuclear receptors
        self.include_nuclear_receptors = include_nuclear_receptors
        # option to print(various parsing statistics)
        self.verbose = verbose
        self.forced = forced
        # inputs dir
        # 2019-08: Updating to use toxcast v3 data
        self.input_dir = "inputs/toxcast-tox21-v3"
        # date present in the file names of the files for this version
        self.version_date = "190708"

        # input files
        self.assay_file = "%s/Assay_Summary_%s.csv"% (self.input_dir, self.version_date) 
        self.s2_assay_file = "%s/S2-ToxCast-Assays.tsv"% (self.input_dir) 
        self.chemical_summary_file = "%s/Chemical_Summary_%s.csv" % (self.input_dir, self.version_date)
        self.zscore_file = "%s/zscore_Matrix_%s.csv" % (self.input_dir, self.version_date)
        self.hitc_file = "%s/hitc_Matrix_%s.csv" % (self.input_dir, self.version_date)
        # instead of using the hitc file, define hits using the ac50 file
        self.ac50_file = "%s/ac50_Matrix_%s.csv" % (self.input_dir, self.version_date)
        self.ac50_cutoff = 50
        #self.chemical_types_file = "%s/chemical_types.tsv" % (self.input_dir)

        # output files
        self.parsed_dir = "%s/parsed" % (self.input_dir)
        t_utils.checkDir(self.parsed_dir)
        self.chem_rec_tfs_file = "%s/chem_rec_tfs.gmt" % (self.parsed_dir) 
        self.chem_hits_file = "%s/chem_prot_hits.csv" % (self.parsed_dir) 
        self.chem_zscore_file = "%s/chem_prot_zscores.csv" % (self.parsed_dir) 

        # mapping dictionaries
        self.chemIDtoName = {}
        self.chemNametoID = {}
        self.chemIDtoTYPE = {}
        self.chemTYPEtoID = defaultdict(set)
        self.assayNametoAccHuman = {}
        self.assayAcctoNameHuman = defaultdict(set)
        self.assayNametoAccListHuman = {}
        self.assayNametoType = {}

        # assay type anaylsis
        self.assay_types = []
        # key is the assay type.
        # each assay type has a list of the proteins perturbed for each chemical 
        self.assay_type_hits = {}
        # key is chemical ID, assay results list (0,1,-1 or NA) is the value
        self.chemical_assay_hits = {}
        # key is chemical ID, list of z-scores is the value
        self.chemical_assay_zscores = {}
        # list of assays in the hits file
        self.hit_assays = []
        # these next two are from the Assay_Summary or S2_Assay_Summary
        # key is assay, value is type_sub
        self.intended_target_type_sub = {} 
        # key is assay, value is family
        self.intended_target_family = {}
        # receptor and tf assays
        # key is assay, value is acc
        self.receptor_assays = {}
        self.tf_assays = {}

        # the sets of hit receptors and TFs (from the Assay summary file)for each chemical
        self.chemical_rec = defaultdict(set)
        self.chemical_tfs = defaultdict(set)
        # proteins stored as uniprot accession IDs
        # prot: 0 or 1. Each protein is labelled a 'hit' if any of the assays are 'hit'
        self.chemical_protein_hit = defaultdict(dict)
        # prot: zscore. This is the largest zscore value of any of the hit (1) assays
        self.chemical_protein_zscore = defaultdict(dict)

        # background interactome networkx object
        #self.net = []

    def get_chemical_maps(self):
        """
        Map the chemical ID to name and other namespaces
        """
        print("reading %s" % (self.chemical_summary_file))
        df = pd.read_csv(self.chemical_summary_file, header=0)
        cols = df.columns
        self.chemIDtoName  = dict(zip(df[cols[3]], df[cols[1]]))
        self.chemNametoID = dict(zip(df[cols[1]], df[cols[3]]))
        #self.chemIDtoDSS   = dict(zip(df[cols[3]], df[cols[0]]))
        #self.chemDSStoID   = dict(zip(df[cols[0]], df[cols[3]]))
        #self.chemDSStoName = dict(zip(df[cols[0]], df[cols[1]]))
        #self.chemNametoDSS = dict(zip(df[cols[1]], df[cols[0]]))

    def load_data(self):
        """
        Load the hits per chemical, zscores of those hits, and the hit rec & tfs per chemical
        """
        self.get_chemical_maps()
        # TODO also write the prots per assay mapping
        if not self.forced and \
                os.path.isfile(self.chem_rec_tfs_file) and \
                os.path.isfile(self.chem_hits_file) and \
                os.path.isfile(self.chem_zscore_file):
            # read directly from the parsed files
            print("reading %s. Use --forced to overwrite" % (self.chem_rec_tfs_file))
            with open(self.chem_rec_tfs_file, 'r') as f:
                for line in f:
                    if line[0] == "#":
                        continue
                    line = line.rstrip().split('\t')
                    chem, prot_type = line[:2]
                    if prot_type == "rec":
                        self.chemical_rec[chem] = line[2:]
                    elif prot_type == "tfs":
                        self.chemical_tfs[chem] = line[2:]
                    else:
                        pass
            print("reading %s. Use --forced to overwrite" % (self.chem_hits_file))
            hits_df = pd.read_csv(self.chem_hits_file, header=0, index_col=0)
            self.chemical_protein_hit = to_dict_drop_na(hits_df.T)
            print("reading %s. Use --forced to overwrite" % (self.chem_zscore_file))
            zscores_df = pd.read_csv(self.chem_zscore_file, header=0, index_col=0)
            self.chemical_protein_zscore = to_dict_drop_na(zscores_df.T)
            # for some reason not all of the hit proteins have a zscore
            # set the default value of the zscore to 0
            for chem, prots in self.chemical_protein_hit.items():
                if chem not in self.chemical_protein_zscore:
                    # defaultdict(int) returns 0 by default
                    self.chemical_protein_zscore[chem] = defaultdict(int)
                    continue
                for p in prots:
                    if p not in self.chemical_protein_zscore[chem]:
                        self.chemical_protein_zscore[chem][p] = 0 

        else:
            # parse all the files
            self.parse_files()

            # now write the output files
            self.write_parsed_files()

    def parse_files(self):
        """
        Parse the chemical summary, assay summary, hit matrix, and zscore matrix, and write them to files
        """
        self.get_chemical_maps()

        # parse the raw files and write the relevant content to file
        print("reading %s" % (self.ac50_file))
        df = pd.read_csv(self.ac50_file, header=0, index_col=0).replace("NA",np.nan).astype(float)
        # header line contains all of the assay names
        self.hit_assays = list(df.columns)
        # keep only the chemicals that have more than 500 non-NA values
        df.dropna(thresh=500, inplace=True)
        print("\t%s chemicals with > 500 non-NAs" % (len(df.index)))
        # define hits to be AC50 < cutoff e.g., 100
        #for ac50_cutoff in [50, 100, 200, 500, 10000]:
        #    print(f"# assays that pass ac50_cutoff {ac50_cutoff}")
        #    print(df.apply(lambda x: x < ac50_cutoff).astype(int).sum().sum())
        print(f"Applying an AC50 cutoff '{self.ac50_cutoff}' to define the hit matrix")
        df_hits = df.applymap(lambda x: x < self.ac50_cutoff if not pd.isnull(x) else x).astype(float)
        self.chemical_assay_hits = df_hits

        print("reading %s" % (self.zscore_file))
        df = pd.read_csv(self.zscore_file, header=0, index_col=0).replace("NA",np.nan)
        print(df.head())
        zscore_assays = list(df.columns)
        # they should match up exactly
        for i in range(len(zscore_assays)):
            if zscore_assays[i] != self.hit_assays[i]:
                print("ERROR: the hit and zscore matrices don't match exactly: %s != %s" % (self.hit_assays[i], zscore_assays[i]))
                sys.exit()
        self.chemical_assay_zscores = df

        print("reading %s" % (self.assay_file))
        # read this one first so that it will be updated with the more recent assay file
        # if anything changed
        # get the assays that are TFs and Receptors from the S2 file received from Richard
        # All but 6 assays of the S2 0 column assay names match the hits file assay name. Using different S2 columns for the assay names didn't make a difference.
        self.read_assay_file(self.s2_assay_file, sep='\t')
        self.read_assay_file(self.assay_file, sep=',')

        # for some reason, most of the 'up' assays have accession numbers, but not the 'down' asays.
        # for example: TG_MRE_CIS_up: Q14872; TG_MRE_CIS_dn: '' 
        # I assume the up and down assays are of the same protein, so I'm fixing that here.
        for assay in self.hit_assays:
            if assay not in self.assayNametoAccHuman and assay[-2:] == 'dn':
                up_assay = assay[:-2] + 'up'
                if up_assay in self.assayNametoAccHuman:
                    #print('correcting %s to %s'%(assay, up_assay))
                    self.assayNametoAccHuman[assay] = self.assayNametoAccHuman[up_assay]
                # also add the dn assay to the receptor/tf dictionaries
                if up_assay in self.receptor_assays and assay not in self.receptor_assays:
                    self.receptor_assays[assay] = self.assayNametoAccHuman[up_assay]
                if up_assay in self.tf_assays and assay not in self.tf_assays:
                    self.tf_assays[assay] = self.assayNametoAccHuman[up_assay]

        # TODO decide when to print(this)
        if self.verbose:
            self.print_assay_summary()
        # account for the assays with multiple accession numbers here
        # for example: ATG_TCF_b_cat_CIS_dn maps to Q9UJU2|P36402|Q9HCS4|Q9NQB0
        for assay in self.assayNametoAccHuman:
            self.assayNametoAccListHuman[assay] = [p for p in self.assayNametoAccHuman[assay].split('|') if p != "NA"]
        # then fix the assayAcctoNameHuman
        assayAcctoNameHuman_fixed = defaultdict(set)
        for prot in self.assayAcctoNameHuman:
            for p in prot.split('|'):
                if p == "NA":
                    continue
                assayAcctoNameHuman_fixed[p].update(self.assayAcctoNameHuman[prot])
        self.assayAcctoNameHuman = assayAcctoNameHuman_fixed

        # don't read in the network here. Write general files first that are non-interactome specific
        #print("reading %s" % (self.interactome))
        #self.net = self.build_network(self.interactome)
        print("building chemical hit receptors, TFs, and proteins")
        self.build_chemical_acc_map()
        #self.get_chemical_types()

    def read_assay_file(self, assay_file, sep=','):
        # now find the receptor and TF assays we will use
        df = pd.read_csv(self.assay_file).replace("NA",np.nan)
        # limit to human assays
        df = df[df['organism'] == "human"]
        # use the uniprot ID from this spreadsheet
        missing_assays = df[~df['assay_component_endpoint_name'].isin(self.hit_assays)]
        if len(missing_assays) > 0:
            print("\t%d missing assays" % (len(missing_assays)))
        df = df[df['assay_component_endpoint_name'].isin(self.hit_assays)]
        self.intended_target_family.update(dict(zip(df['assay_component_endpoint_name'], df['intended_target_family'])))
        self.intended_target_type_sub.update(dict(zip(df['assay_component_endpoint_name'], df['intended_target_type_sub'])))
        # only keep the assays which have a uniprot ID
        # TODO map the gene names to uniprot IDs for the other assays?
        #df = df[df['intended_target_uniprot_accession_number'] != np.nan]
        df.dropna(axis=0, subset=['intended_target_uniprot_accession_number'], inplace=True)
        # remove bioseek assays
        df = df[~df['assay_component_endpoint_name'].str.contains("BSK")]
        df_rec = df[df['intended_target_type_sub'] == 'receptor']
        df_tfs = df[df['intended_target_type_sub'] == 'transcription factor']
        if not self.include_nuclear_receptors:
            df_rec = df_rec[df_rec['intended_target_family'] != "nuclear receptor"]
            df_tfs = df_tfs[~df_tfs['intended_target_family'].isin(["nuclear receptor", "gpcr"])]
        self.receptor_assays.update(dict(zip(df_rec['assay_component_endpoint_name'], 
            df_rec['intended_target_uniprot_accession_number'])))
        # remove attagene assays as those are TF assays (none of them are hit for receptors anyway)
        self.receptor_assays = {a: prots for a, prots in self.receptor_assays.items() if 'ATG' not in a}
        self.tf_assays.update(dict(zip(df_tfs['assay_component_endpoint_name'], 
            df_tfs['intended_target_uniprot_accession_number'])))
        #self.receptors = set(p for receptors in self.receptor_assays.values() for p in receptors)
        #self.tfs = set(p for tfs in self.tf_assays.values() for p in tfs)
        # also store the prots for each assay, and the assays for each prot
        self.assayNametoAccHuman.update(dict(zip(df['assay_component_endpoint_name'], 
            df['intended_target_uniprot_accession_number'])))
        for acc, assay_name in zip(df['intended_target_uniprot_accession_number'], 
                                   df['assay_component_endpoint_name']):
            self.assayAcctoNameHuman[acc].add(assay_name)

    def print_assay_summary(self):
        print("%s assays"%str(len(self.hit_assays)))
    #    print("%s human assays"%str(len(human_assays)))
    #    print("%s s2 human assays"%str(len(s2_human_assays)))
    #    print("%s human assays & s2 human assays"%str(len(human_assays & s2_human_assays)))
    #    print("%s human assays - s2 human assays"%str(len(human_assays - s2_human_assays)))
    #    print("%s s2 human assays - human assays"%str(len(s2_human_assays - human_assays)))
        print("%s human assays with accession numbers"%str(sum([1 for assay in self.hit_assays if assay in self.assayNametoAccHuman])))
        print("%s human accession numbers"%str(len(self.assayAcctoNameHuman)))
        print("%s receptor assays"%str(len(self.receptor_assays)))
        print("%s tf assays"%str(len(self.tf_assays)))
        print("%s receptor acc"%str(len(set(p for receptors in self.receptor_assays.values() for p in receptors))))
        print("%s tf acc"%str(len(set(p for tfs in self.tf_assays.values() for p in tfs))))
        # skip this for now as it doesn't seem to be working
        #print("Checking to see if all acc are primary")
        #csbdb_interface = csbdb.Interface()
        #for accs, assays in self.assayAcctoNameHuman.items():
        #    #print(accs, assays)
        #    if len(accs.split('|')) > 1:
        #        print('multiple accs: ', accs)
        #    for acc in accs.split('|'):
        ##        primary_acc = csbdb_interface.make_primary(acc, taxonomy_number='9606')
        #        if len(primary_acc) != 1 or list(primary_acc)[0] != acc:
        #            print("%s not primary. mapped to %s" % (acc, str(primary_acc)))

    def build_chemical_acc_map(self, zscore_cutoff=None):
        """ zscore_cutoff is an option to only include the hits that have a zscore > X (e.g., 3)
        """
        chemical_protein_hit_list = {}
        # structure of the dictionary being built:
        # keep track of the uniprot IDs not in the network
        #not_in_net = set()
        for i, chemical in enumerate(tqdm(self.chemical_assay_hits.index)):
            curr_chem_prot_hit_list = defaultdict(list)
            curr_chem_prot_hit = defaultdict(list)
            # we are looping through the assays, but what we really want to determine is if the 
            # proteins are responsive to the chemical, or hit.
            # as such, we find all of the assays corresponding to a protein 
            for assay in self.hit_assays:
                #print(i, len(self.chemical_assay_hits[chemical]))
                hit = self.chemical_assay_hits[assay].iloc[i]
                zscore = self.chemical_assay_zscores[assay].iloc[i]
                #zscore = float(zscore) if zscore != "NA" else zscore
                # if this assay was successfully tested (0 or 1) and the assay is of a human protein
                if hit != "NA" and not np.isnan(hit) and \
                        int(hit) >= 0 and assay in self.assayNametoAccListHuman: 
                    prots = self.assayNametoAccListHuman[assay]
                    for prot in prots:
                        # Many proteins have both an up/down version, so if either of those is hit,
                        # we will set it as a hit for the chemical
                        # 2023-02-02 UPDATE: For the 7 proteins tested by Tox21 assays, they overlap with ATG assays.
                        # Only label a protein as hit if there is consensus among all the assays testing that protein.
                        # To determine a consensus, append the hit values to a list
                        curr_chem_prot_hit_list[prot].append((assay, int(hit), assay in self.receptor_assays, assay in self.tf_assays))
                        if int(hit) == 1 and (zscore_cutoff is None or zscore > zscore_cutoff):
                            # keep the biggest zscore for this protein
                            if prot not in self.chemical_protein_zscore[chemical] or zscore > self.chemical_protein_zscore[chemical][prot]:
                                self.chemical_protein_zscore[chemical][prot] = zscore
                            curr_chem_prot_hit[prot].append((1, 'TOX21' in assay))
                            # add the assay's corresponding protein to the list of receptor or tf assays
                            # don't allow duplicates
                            if assay in self.receptor_assays:
                                self.chemical_rec[chemical].add(prot)
                            if assay in self.tf_assays: 
                                self.chemical_tfs[chemical].add(prot)
                        else: 
                            curr_chem_prot_hit[prot].append((0, 'TOX21' in assay))
            # for each protein, determine if it is a hit or not
            chem_prot_hit = {}
            for prot, hits_and_types in curr_chem_prot_hit.items():
                non_tox21_hits = [hit for hit, tox21 in hits_and_types if not tox21]
                tox21_hits = [hit for hit, tox21 in hits_and_types if tox21]
                if len(non_tox21_hits) > 0 and len(tox21_hits) > 0:
                    hit = min([max(non_tox21_hits)] + tox21_hits)
                elif len(non_tox21_hits) > 0:
                    hit = max(non_tox21_hits)
                elif len(tox21_hits) > 0:
                    hit = min(tox21_hits)
                chem_prot_hit[prot] = hit
            chemical_protein_hit_list[chemical] = curr_chem_prot_hit_list
            self.chemical_protein_hit[chemical] = chem_prot_hit
            # update the chemical's receptors and tfs
            self.chemical_rec[chemical] = set(p for p in self.chemical_rec[chemical] if chem_prot_hit[p] == 1)
            self.chemical_tfs[chemical] = set(p for p in self.chemical_tfs[chemical] if chem_prot_hit[p] == 1)

        # now write the hits to a file to see how often a protein is only a "hit" for a small fraction of assays
        out_file = "%s/chem_prot_assay_hit_vals.csv" % (self.parsed_dir) 
        with open(out_file, 'w') as out:
            out.write("chem,prot,assay,hit\n")
            for chem, prot_hit_list in chemical_protein_hit_list.items():
                for prot, assay_hit_list in prot_hit_list.items():
                    for assay, hit, rec_assay, tf_assay in assay_hit_list:
                        out.write(f"{chem},{prot},{assay},{hit},{rec_assay},{tf_assay}\n")

        #print("\t%s uniprot IDs were not in the network and were ignored: %s" % (len(not_in_net), ', '.join(not_in_net)))
        # don't allow proteins to be in both the list of receptors and the list of tfs
        if not self.include_nuclear_receptors:
            removed = set()
            for chemical in self.chemical_rec:
                new_rec = self.chemical_rec[chemical] - self.chemical_tfs[chemical]
                new_tfs = self.chemical_tfs[chemical] - self.chemical_rec[chemical]
                removed.update(self.chemical_rec[chemical] & self.chemical_tfs[chemical])
                self.chemical_rec[chemical] = new_rec
                self.chemical_tfs[chemical] = new_tfs
            print("\t%s prot(s) listed as both receptors and TFs. Removing them from analysis since they are likely nuclear receptors." % (','.join(removed)))

    # get the type of each chemical (i.e., Bactericide, Pharmaceutical, Solvent, etc)
    def get_chemical_types(self):
        print("reading %s" % (self.chemical_types_file))
        with open(self.chemical_types_file, 'r') as tsv:
            for line in tsv:
                if line[0] == '#':
                    continue
                line = line.rstrip().split('\t')
                chemID = line[0]
                chemical_type = line[2].replace(' ', '_')
                self.chemIDtoTYPE[chemID] = chemical_type
                self.chemTYPEtoID[chemical_type].add(chemID)

    def write_chem_rec_tfs(self, out_file):
        print("writing %s" % (out_file))
        with open(out_file, 'w') as out:
            out.write("#chem\ttype\tprots\n")
            for chem in self.chemical_rec:
                rec = sorted(self.chemical_rec[chem])
                tfs = sorted(self.chemical_tfs[chem])
                out.write("%s\trec\t%s\n" % (chem, '\t'.join(rec)))
                out.write("%s\ttfs\t%s\n" % (chem, '\t'.join(tfs)))

    def write_parsed_files(self):
        self.write_chem_rec_tfs(self.chem_rec_tfs_file)

        # now build a dataframe of the hits and zscore
        print("writing %s" % (self.chem_hits_file))
        hits_df = pd.DataFrame(self.chemical_protein_hit)
        # sort by protein ID(?)
        hits_df = hits_df.sort_index().T
        hits_df.to_csv(self.chem_hits_file)
        print("writing %s" % (self.chem_zscore_file))
        zscores_df = pd.DataFrame(self.chemical_protein_zscore)
        # sort by protein ID(?)
        zscores_df = zscores_df.sort_index().T
        zscores_df.to_csv(self.chem_zscore_file)

    def limit_to_interactome(self, G):
        """
        Limit the hit receptors and TFs to these nodes. 
        Also limit the chemicals to those that have at least 1 hit receptor and TF.
        """
        nodes = set(G.nodes())
        # minimum number of hit rec and tfs for each chemical
        min_num = 1
        # keep track of the receptors with no outgoing edges
        rec_no_out = set()
        # and the tfs with no incoming edges
        tfs_no_in = set()
        orig_num_chem = len(self.chemical_rec)
        new_chem_rec, new_chem_tfs = {}, {}
        for c, rec in self.chemical_rec.items():
            rec_in_net = set(rec) & nodes
            for n in rec_in_net.copy():
                if G.out_degree(n) == 0:
                    rec_in_net.discard(n)
                    rec_no_out.add(n)
            if len(rec_in_net) == 0:
                continue
            new_chem_rec[c] = rec_in_net
        for c, tfs in self.chemical_tfs.items():
            tfs_in_net = set(tfs) & nodes
            for n in tfs_in_net.copy():
                if G.in_degree(n) == 0:
                    tfs_in_net.discard(n)
                    tfs_no_in.add(n)
            if len(tfs_in_net) == 0:
                continue
            new_chem_tfs[c] = tfs_in_net
        # now limit to the chemicals with at least 1 hit rec and 1 hit tf
        common_chems = set(new_chem_rec.keys()) & set(new_chem_tfs.keys())
        self.chemical_rec = {}
        self.chemical_tfs = {}
        for chem in common_chems:
            rec, tfs = new_chem_rec[chem], new_chem_tfs[chem]
            if len(rec) < min_num or len(tfs) < min_num:
                continue
            self.chemical_rec[chem] = rec
            self.chemical_tfs[chem] = tfs
        print("\t%d chemicals with at least 1 hit receptor & TF" % (len(self.chemical_rec)))
        print("\t%d total rec, %d total tfs."  % (
            len(set(r for rec in self.chemical_rec.values() for r in rec)),
            len(set(t for tfs in self.chemical_tfs.values() for t in tfs))))
        print("\t%d unreachable rec, %d unreachable tfs" % (
            len(rec_no_out), len(tfs_no_in)))

    def write_num_chemicals_per_hit(self, out_file):
        # dictionary of protein hits: set of chemicals
        protein_hits_chemicals = {}
        for chemical in self.chemicals:
            if len(self.chemicals[chemical]['rec']) >= 1 and len(self.chemicals[chemical]['tfs']) >= 1:
                protein_hits_chemicals[chemical] = self.chemical_protein_hit[chemical]
        # write the number of chemicals that perturb each protein
        with open(out_file, 'w') as out:
            out.write('\n'.join(["%s\t%d" % (prot, protein_hits_chemicals[prot]) for prot in protein_hits_chemicals]) + '\n')

    def write_chem_num_rectfs(self, out_file):
        print("Writing # of rec and tfs perturbed by each chemical to %s" % out_file)
        with open(out_file, 'w') as out:
            out.write("#Chemical Name\tID\t# rec\t# tfs\n")
            for chem in self.chemical_rec:
                rec, tfs = self.chemical_rec[chem], self.chemical_tfs[chem]
                out.write('\t'.join([self.chemIDtoName[chem], chem, str(len(rec)), str(len(tfs))]) + '\n')

    def write_chemical_perturbed_rec_tfs(self, chemicals_file, rec_tfs_dir, include_zscore_weight=False):
        """ write the chemicals file as well as the perturbed rec and tfs for each chemical
        We are writing a single file for each chemical so we can use each file for running pathlinker/cyclinker
        """
        print("\tWriting the chemicals file (%s) as well as the perturbed rec and tfs for each chemical in %s" % (chemicals_file, rec_tfs_dir))
        t_utils.checkDir(rec_tfs_dir)

        # first write a general file with the hit rec and tf per chemical
        out_file = "%s/../chem_rec_tfs.gmt" % (rec_tfs_dir)
        self.write_chem_rec_tfs(out_file)
        # also write a table with the number of hit rec and tfs per chemical
        out_file = "%s/../chemical_num_rectfs.txt" % (rec_tfs_dir)
        self.write_chem_num_rectfs(out_file)

        # keep track of all of the receptors and tfs and write them to a file as well
        all_rec = set()
        all_tfs = set()

        chemicals = self.chemical_rec.keys()

        # first write the chemicals
        with open(chemicals_file, 'w') as out:
            out.write('\n'.join(["%s\t%s" % (chemical, self.chemIDtoName[chemical]) for chemical in chemicals]))

        if include_zscore_weight is True:
            zscores = []
            for chem, prots in self.chemical_rec.items():
                zscores += [self.chemical_protein_zscore[chem][p] for p in prots]
            for chem, prots in self.chemical_tfs.items():
                zscores += [self.chemical_protein_zscore[chem][p] for p in prots]
            zscores = [z for z in zscores if not np.isnan(z)]
            # use the maximum to normalize the zscores
            max_zscore = max(zscores)
            print("max zscore is: %0.2f" % max_zscore)

        for chem in tqdm(chemicals):
            # some of the self.chemicals have spaces in their names, so use the ID rather than the name.
            rec = set(self.chemical_rec[chem])
            tfs = set(self.chemical_tfs[chem])
            all_rec.update(rec)
            all_tfs.update(tfs)
            if include_zscore_weight is False:
                t_utils.writeRecTFs("%s/%s-rec-tfs.txt" % (rec_tfs_dir, chemical), rec, tfs)
            else:
                # convert the zscore to a cost by taking 1 - (zscore / max zscore)
                # the lower the zscore, the higher the cost will be
                zscores = {}
                curr_zscores = self.chemical_protein_zscore[chem]
                for prots in (self.chemical_rec[chem], self.chemical_tfs[chem]):
                    for p in prots:
                        zscore = curr_zscores[p] 
                        zscores[p] = zscore if not pd.isnull(zscore) else 0
                costs = {p: 1-(zscore / float(max_zscore)) for p, zscore in zscores.items()}
                t_utils.writeRecTFs("%s/%s-rec-tfs.txt" % (rec_tfs_dir, chem), rec, tfs, costs=costs, zscores=zscores)

        out_file = "%s/all-rec-tfs.txt" % (rec_tfs_dir)
        print("Writing all of the assayed receptors and tfs to the file: %s" % (out_file))
        t_utils.writeRecTFs(out_file, all_rec, all_tfs)

    # @param json_data the json dictionary to be written
    # @param path the path of the json file
    def write_json(self, json_data, file_path):
        with open(file_path, 'w') as newJsonFile:
            json.dump(json_data, newJsonFile, sort_keys=True, indent=4)

    # How many assay accession numbers are in the human receptors and tfs list
    def map_assays_rec_tfs(self):
        #for chemical, assays in self.chemical_assay_hits.iteritems():
        #    break

        receptors = set()
        with open("pathlinker-data/human-receptors-acc.txt", 'r') as rec_in:
            for line in rec_in:
                if line[0] != "#":
                    receptors.add(line.strip())
        tfs = set()
        with open("pathlinker-data/human-tfs-acc.txt", 'r') as tfs_in:
            for line in tfs_in:
                if line[0] != "#":
                    tfs.add(line.strip())

        toxcast_receptors = set()
        toxcast_tfs = set()
        # We souldn't need to map the assay names
        for assay in self.hit_assays:
            if assay in self.assayNametoAccHuman:
                if self.assayNametoAccHuman[assay] in tfs:
                    #print("%s - %s is a human transcription factor"%(assay, self.assayNametoAccHuman[assay]))
                    toxcast_tfs.add(self.assayNametoAccHuman[assay])
                elif self.assayNametoAccHuman[assay] in receptors:
                    #print("%s - %s is a human transcription factor"%(toxcast_acc, name))
                    toxcast_receptors.add(self.assayNametoAccHuman[assay])
        print("%s assays map to human receptors"%(len(toxcast_receptors)))
        print("%s assays map to human tfs"%(len(toxcast_tfs)))
        self.build_network(self.interactome,toxcast_receptors,toxcast_tfs)

    def write_assayed_rec_tfs(self, rec_tfs_files, out_file):
        """ write the complete list of receptors and tfs from the intput file that are in the interactome
        Should only be called after each chemical's list of rec and tfs have been written
        """
        print("\tCreating the '%s' file set of receptors and tfs from %d rec_tfs_files" % (out_file, len(rec_tfs_files)))
        if not self.net:
            self.net = self.build_network(self.interactome)

        receptors = set()
        tfs = set()
        for rec_tfs_file in rec_tfs_files:
            # first get all of the receptors and tfs
            rec, t = t_utils.getRecTFs(rec_tfs_file)
            receptors.update(set(rec))
            tfs.update(set(t))

        # all of the rec and tfs should be in the interactome
        # No protein should be both a rec and tf
        # remove receptors and tfs that aren't in the interactome
        orig_len_rec = len(receptors)
        orig_len_tfs = len(tfs)
        receptors = set([r for r in receptors if r in self.net and len(self.net.out_edges(r)) > 0])
        tfs = set([tf for tf in tfs if tf in self.net and len(self.net.in_edges(tf)) > 0])
        print("\tRemoved %d recptors and %d tfs not in the interactome" % (orig_len_rec - len(receptors), orig_len_tfs - len(tfs)))
        # remove receptors and tfs that are in both
        orig_len_rec = len(receptors)
        orig_len_tfs = len(tfs)
        receptors = set([r for r in receptors if r not in tfs])
        tfs = set([tf for tf in tfs if tf not in receptors])
        print("\tRemoved %d recptors and %d tfs that were in both the receptors and tfs sets" % (orig_len_rec - len(receptors), orig_len_tfs - len(tfs)))

        # now write the output file
        t_utils.checkDir(os.path.dirname(out_file))
        t_utils.writeRecTFs(out_file, receptors, tfs)

    def write_global_rec_tfs(self, rec_tfs_file, out_file):
        """ write the complete list of receptors and tfs from the intput file that are in the interactome
        """
        print("\t Creating the '%s' file set of receptors and tfs from %s" % (out_file, rec_tfs_file))
        if not self.net:
            self.net = self.build_network(self.interactome)

        # first get the list of all human rec and tfs
        receptors, tfs = t_utils.getRecTFs(rec_tfs_file)

        # all of the rec and tfs should be in the interactome AND rec have outgoing edges, tfs incoming edges
        # No protein should be both a rec and tf
        # remove receptors and tfs that aren't in the interactome
        receptors = set([r for r in receptors if r in self.net and len(self.net.out_edges(r)) > 0])
        tfs = set([tf for tf in tfs if tf in self.net and len(self.net.in_edges(tf)) > 0])
        # remove receptors and tfs that are in both
        receptors = set([r for r in receptors if r not in tfs])
        tfs = set([tf for tf in tfs if tf not in receptors])

        # now write the output file
        t_utils.checkDir(os.path.dirname(out_file))
        t_utils.writeRecTFs(out_file, receptors, tfs)

    def write_all_rec_tfs(self, out_file):
        """ write the complete list of receptors and tfs perturbed by any chemical as well as the assay name
        """
        print("\t Creating the '%s' file set of receptors and tfs" % (out_file))
        if not self.net:
            self.net = self.build_network(self.interactome)

        t_utils.checkDir(os.path.dirname(out_file))
        out = open(out_file, 'w')
        out.write("#uniprot_acc\tnode_type\tassay_name\n")
        for r in self.receptor_assays:
            for acc in self.assayNametoAccListHuman[r]:
                if acc in self.net:
                    out.write(acc+'\t'+'receptor'+'\t'+r+'\n')
        for tf in self.tf_assays:
            for acc in self.assayNametoAccListHuman[tf]:
                if acc in self.net:
                    out.write(acc+'\t'+'tf'+'\t'+tf+'\n')
        out.close()

#def write_chemicals_to_file(out_file):
#    out = open(out_file, 'w')
#    ids = []
#    with open("toxcast-data/2013_Data.csv") as data:
#        # skip the header lines
#        for i in range(0,4): data.readline()
#        # line containing the self.chemicals
#        self.chemicals_header = data.readline().strip().split('\t')
#        data.readline()
#        ids = data.readline().strip().split('\t')[2:]
#    count = 0
#    for chemical in sorted(self.chemical_assay_hits):
#        chemical_match = '%s\t%s\t'%(self.chemIDtoName[chemical],self.chemIDtoDSS[chemical])
#        for DSSID in ids:
#            id_num = DSSID.split("_")[-1]
#            if self.chemIDtoDSS[chemical] == id_num:
#                count += 1
#                chemical_match += id_num
#        out.write(chemical_match+'\n')
#    print("%s self.chemicals matched"%(str(count)))
#    out.close()


def to_dict_drop_na(df):
    # code to drop na values when converting to dictionary
    # inspired by this stack overflow answer: https://stackoverflow.com/a/46098323
    df_dict = df.to_dict()
    df_dict = {row: {col:v for col,v in col_vals.items() if pd.notnull(v)} for row, col_vals in df_dict.items()}
    return df_dict


def parse_args(args):
    parser = OptionParser()
    #parser.add_option('-t','--threshold', type='float', metavar='STR', default=0, help='the assay result threshold to consider (default = 0)')
    #parser.add_option('-o','--output', type='string', metavar='STR', default='chemical-rec-tfs.txt', help='the output file to write to (default = 0)')
    #parser.add_option('-O','--output_dir', type='string', metavar='STR', default='results/self.chemicals/pathlinker-input/', help='the output dir to write pathlinker input to')
    #parser.add_option('','--num_rec', type='int', metavar='STR', default=1, help='min # of receptors perturbed by a chemical (default = 1)')
    #parser.add_option('','--num_tfs', type='int', metavar='STR', default=0, help='min # of tfs perturbed by a chemical. If 0, will write all human tfs as targets. (default = 0)')
    #parser.add_option('','--interactome',action='store',  # default='pathlinker-data/2015pathlinker-weighted.txt',\
    #        help='Background interactome used when running pathlinker')
    #parser.add_option('','--write-hit-proteins', type='string',\
    #        help='Write each chemicals perturbed proteins to a specified dir')
    #parser.add_option('','--write-num-chemicals-per-hit', type='string',\
    #        help='Write the number of chemicals that hit each protein')
    parser.add_option('','--forced', action='store_true', default=False,
            help='Force reading the data files from scratch and overwriting the parsed files')
    #parser.add_option('','--all-tfs', type='store_true', metavar='STR', default=1, help='run pathlinker using all human tfs')
    #parser.add_option('-a','--no_activator', type='float', metavar='STR', default=0, help='the assay result threshold to consider (default = 0)')

    (opts, args) = parser.parse_args()
    #print(opts)

    return opts


if __name__ == '__main__':
    # TODO add an interactome argument?
    opts = parse_args(sys.argv)
    print(opts.forced)
    #ToxCastData(verbose=True, forced=opts.forced).main(opts.write_hit_proteins, opts.write_num_chemicals_per_hit)
    ToxCastData(verbose=True, forced=opts.forced).load_data()
    #main()

