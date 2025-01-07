"""
Helper functions for analyzing fastsimbac output.
"""

import numpy as np
import pandas as pd

class fastsimbacHelper:
    def __init__(self, fsb_path, sim_params):
        self.fsb_path = fsb_path
        self.sim_params = sim_params
        self.snp_df = self.load_fsb_data(fsb_path)

    @staticmethod
    def load_fsb_data(filename):
        """
        Load fastsimbac output haplotype data.

        Columns are:
        - pos: position of the polymorphic site (float between 0 and 1)
        - freq: frequency of the derived allele
        - haplotype columns for each individual (0 or 1)
        """
        f = open(filename)
        all_sites = []
        haplotypes = []
        for line in f:
            if line.startswith('SITE'):
                items = line.split()
                all_sites.append(items[2:4])
                haplotypes.append(items[4])
        f.close()
        # turn the string of 0s and 1s into a list of integers
        haplotypes = [list(map(int, list(h))) for h in haplotypes]
        snp_info = pd.DataFrame(all_sites, columns=['pos', 'freq'], 
                                dtype=float)
        haplotypes = pd.DataFrame(haplotypes)

        snp_df = pd.concat([snp_info, haplotypes], axis=1)
        return snp_df
    
    def get_pair_snps(self, pair):
        """
        Get the SNPs that differ between two individuals in a pair.
        Return a slice of the snp_df with only the relevant columns.
        """
        ind1, ind2 = pair
        sub_df = self.snp_df.loc[(self.snp_df[ind1] != self.snp_df[ind2])]
        return sub_df[['pos', ind1, ind2]].copy()
    
    def get_pair_identical_runs(self, pair):
        """
        Get the runs of identical SNPs between two individuals in a pair.

        Returns a list of run lengths
        Note this is in fraction of genome length, not number of sites
        """
        snps = self.get_pair_snps(pair)
        # pad two snps at the beginning and end to catch runs at the ends
        new_row_top = pd.DataFrame([[0, 0, 1]], columns=['pos', pair[0], pair[1]])
        new_row_bottom = pd.DataFrame([[1, 1, 0]], columns=['pos', pair[0], pair[1]])
        snps = pd.concat([new_row_top, snps, new_row_bottom], ignore_index=True)
        runs = snps['pos'].diff().dropna().values
        return runs
    
    def get_max_runs(self, pairs):
        max_runs = []
        for pair in pairs:
            runs = self.get_pair_identical_runs(pair)
            max_runs.append(runs.max())
        return max_runs

"""
Helper functions for generating pairs of individuals for comparison.
""" 

def generate_pairs(samples):
    # helper for doing pairwise comparisons
    # generate all unique pairs within a list
    pairs = []
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            pairs.append((samples[i], samples[j]))
    return pairs

def generate_pairs_between(samples1, samples2):
    # generate all unique pairs between two lists
    pairs = []
    for i in range(len(samples1)):
        for j in range(len(samples2)):
            pairs.append((samples1[i], samples2[j]))
    return pairs

def sample_pairs(pair_set, n):
    # sample n pairs from a set of pairs
    sampled = pair_set[np.random.choice(np.arange(len(pair_set)), size=n, replace=False)]
    return sampled

def sample_within_between_pairs(pop_size, sample_size):
    pop1 = np.arange(0, pop_size)
    pop2 = np.arange(pop_size, 2*pop_size)
    pop1_pairs = generate_pairs(pop1)
    pop2_pairs = generate_pairs(pop2)
    all_pairs = pop1_pairs + pop2_pairs
    indices = np.random.choice(len(all_pairs), size=sample_size, replace=False)
    sampled_within_pairs = [all_pairs[i] for i in indices]

    # sample pairs between sets of individuals
    pop12_pairs = generate_pairs_between(pop1, pop2)
    indices = np.random.choice(len(pop12_pairs), size=sample_size, replace=False)
    sampled_between_pairs = [pop12_pairs[i] for i in indices]
    return sampled_within_pairs, sampled_between_pairs