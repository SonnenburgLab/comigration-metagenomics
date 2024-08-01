import numpy as np
import pandas as pd
import scipy
from scipy.cluster.hierarchy import linkage, fcluster
from utils import metadata_utils

import config

def load_hgt_res(hgt_res_path=config.hgt_res_path):
    hgt_res = pd.read_csv(hgt_res_path)
    hgt_res['species'] = hgt_res['name'].str.replace('_v7', '')
    hgt_res['genome1'] = hgt_res['querry'].str.replace('\.fa', '')
    hgt_res['genome2'] = hgt_res['reference'].str.replace('\.fa', '')
    return hgt_res

def long_form_to_symmetric(df, row_name='genome1', col_name='genome2', val_name='div'):
    # Use combine_first to fill NaN values in the pivot_table with those in transpose_table
    pivot_table = df.pivot(index=row_name, columns=col_name, values=val_name)
    transpose_table = pivot_table.T
    symmetric_table = pivot_table.combine_first(transpose_table)
    np.fill_diagonal(symmetric_table.values, 0)
    return symmetric_table


def cluster_close_pairs(pd_mat, cut_dist=1e-3, return_clusters=False):
    """
    Take a pairwise divergence matrix (square form, index & column are MAG names)
    and cluster the MAGs based on their pairwise divergence.
    Returns a list of MAG names
    """
    condensed_y = scipy.spatial.distance.squareform(pd_mat)
    Z = linkage(condensed_y, method='average')

    clusters = fcluster(Z, cut_dist, criterion='distance')
    cluster_samples = pd.DataFrame(index=pd_mat.index, data={'cluster': clusters})
    chosen_samples = sample_genome_per_cluster(cluster_samples)
    if return_clusters:
        return chosen_samples, cluster_samples
    return chosen_samples

def sample_genome_per_cluster(cluster_samples):
    # Iterate over each cluster and select one sample
    chosen_samples = []
    for cluster in cluster_samples['cluster'].unique():
        cluster_data = cluster_samples[cluster_samples['cluster'] == cluster]
        # Select one random sample from each cluster
        sample = cluster_data.sample(1)
        chosen_samples.append(sample)
    # Concatenate all selected samples into a single DataFrame
    chosen_samples_df = pd.concat(chosen_samples).sort_index()
    return chosen_samples_df.index

def filter_df_by_samples(df, chosen_samples):
    mask = df['genome1'].isin(chosen_samples) & df['genome2'].isin(chosen_samples)
    return df[mask]

def dedup_pairwise_df(df, cut_dist=1e-3):
    symmetric_table = long_form_to_symmetric(df, row_name='genome1', col_name='genome2', val_name='div')
    dedup_genomes = cluster_close_pairs(symmetric_table, cut_dist=cut_dist)
    dedup_df = filter_df_by_samples(df, dedup_genomes)
    return dedup_df

def compute_pop_comparisons(df, stat_func):
    """
    Compute the mean statistic between populations for each pairwise comparison

    df: DataFrame
        A long-form DataFrame with columns 'pop1', 'pop2', 'div' and more
    stat_func: function
        A function that takes a DataFrame and returns a statistic
        e.g. lambda x: x['div'].mean()
    """
    included_pops = np.unique(df[['pop1', 'pop2']])
    pop_comp_mat = pd.DataFrame(columns=included_pops, index=included_pops, dtype=float)

    for comp, x in df.groupby('comp'):
        stat = stat_func(x)
        pop1, pop2 = comp.split('-')[0], comp.split('-')[1]
        pop_comp_mat.loc[pop1, pop2] = stat
        pop_comp_mat.loc[pop2, pop1] = stat
    return pop_comp_mat


"""
Below are some functions for summarizing pairwise statistics into species-comparison summaries
Current ideas as of 2024-07-14:
- HGT_score: compare quantile of max_run to a reference length (e.g. 700bp)
- long_run_length: simply the quantile length
"""

def HGT_score(df, quantile=0.99, ref_length=700):
    """
    Idea: take the distribution of max_runs, find the tail length at a certain quantile,
    then compare that with a reference length to get a score

    For example, 700bp corresponds to a length scale of 1 SNV in 10,000 years
    """
    hgt_len = df['max_run'].quantile(quantile, interpolation='higher')
    hgt_score = np.log10(hgt_len / ref_length)
    return hgt_score


def long_run_length(df, quantile=0.99):
    return df['max_run'].quantile(quantile, interpolation='higher')


class PairwiseHelper:
    """
    A class to help with pairwise comparisons
    TODO: eventually add the identical gene and genome ANI statistics
    TODO: might need to compute identical block myself because drep cannot finish every pair
    """
    def __init__(self, databatch):
        self.databatch = databatch
        self.metadata = metadata_utils.MetadataHelper(databatch)
        self.pairwise_summary = self.load_pairwise_summary()
        self.hgt_summary = load_hgt_res()

    def load_pairwise_summary(self):
        pairwise_summary_path = config.run_path / f'{self.databatch}_annotated.csv'
        return pd.read_csv(pairwise_summary_path)

    def get_species_pairwise_summary(self, species_name):
        species_res = self.pairwise_summary[self.pairwise_summary['species'] == species_name]
        return species_res

    def filter_close_pairs_by_div(self, species_res):
        average_div = species_res['div'].mean()
        cut_dist = average_div / config.close_pair_div_ratio
        return dedup_pairwise_df(species_res, cut_dist=cut_dist)

    def filter_clonal_pairs_by_perc_id(self, species_res, threshold=config.clonal_cluster_pi_threshold):
        passed_mags = self.get_nonclonal_mags(species_res['species'].iloc[0], threshold)
        return filter_df_by_samples(species_res, passed_mags)
    
    def get_nonclonal_mags(self, species_name, threshold=config.clonal_cluster_pi_threshold):
        """
        Filter the species pairwise comparison by clonal cluster
        Requires dRep HGT results to contain this species

        Any pair above the threshold will be de-duplicated
        """
        species_hgt_res = self.hgt_summary[self.hgt_summary['species'] == species_name]
        symmetric_table = long_form_to_symmetric(species_hgt_res, row_name='genome1', col_name='genome2', val_name='perc_id')
        np.fill_diagonal(symmetric_table.values, 1)
        symmetric_table = 1 - symmetric_table

        passed_mags = cluster_close_pairs(symmetric_table, 1-threshold)
        return passed_mags

    def filter_population(self, species_res, pop):
        return species_res[~((species_res['pop1'] == pop) | (species_res['pop2'] == pop))]
