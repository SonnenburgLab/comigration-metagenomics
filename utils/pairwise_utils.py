import numpy as np
import pandas as pd
import scipy
import pickle
from scipy.cluster.hierarchy import linkage, fcluster
from utils import metadata_utils

import config

def load_hgt_res(hgt_res_path=config.hgt_res_path):
    hgt_res = pd.read_csv(hgt_res_path)
    hgt_res['species'] = hgt_res['name'].str.replace('_v7', '')
    hgt_res['genome1'] = hgt_res['querry'].str.replace('.fa', '')
    hgt_res['genome2'] = hgt_res['reference'].str.replace('.fa', '')
    return hgt_res

def long_form_to_symmetric(df, row_name='genome1', col_name='genome2', val_name='div', diagonal_fill=0):
    """
    Convert a long-form DataFrame to a symmetric DataFrame
    Meant to deal with pairwise comparisons
    """
    # Use combine_first to fill NaN values in the pivot_table with those in transpose_table
    pivot_table = df.pivot(index=row_name, columns=col_name, values=val_name)
    transpose_table = pivot_table.T
    symmetric_table = pivot_table.combine_first(transpose_table)
    np.fill_diagonal(symmetric_table.values, diagonal_fill)
    return symmetric_table

def linkage_clustering(pd_mat):
    """
    Take a pairwise divergence matrix (square form, index & column are MAG names)
    and cluster the MAGs based on their pairwise divergence.
    """
    condensed_y = scipy.spatial.distance.squareform(pd_mat)
    Z = linkage(condensed_y, method='average')
    return Z

def cluster_close_pairs(pd_mat, cut_dist=1e-3):
    """
    Take a pairwise divergence matrix (square form, index & column are MAG names)
    and cluster the MAGs based on their pairwise divergence.
    """
    Z = linkage_clustering(pd_mat)
    clusters = fcluster(Z, cut_dist, criterion='distance')
    cluster_samples = pd.DataFrame(index=pd_mat.index, data={'cluster': clusters})
    return cluster_samples

def sample_genome_per_cluster(cluster_samples):
    """
    Given a DataFrame with columns 'cluster' and index as MAG names,
    select one sample from each cluster and return a list of
    MAG names
    """
    # Iterate over each cluster and select one sample
    # Returns a list of MAG names
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
def compute_L99(df, quantile=0.99):
    l99 = df['max_run'].quantile(quantile, interpolation='higher')
    return l99

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


"""
Below are some classes for organizing pairwise comparisons
"""
class SpeciesPairwiseHelper:
    def __init__(self, species_name, run_summary, hgt_summary, metadata, cluster_threshold=config.clonal_cluster_pi_threshold):
        self.metadata = metadata
        self.species_name = species_name
        self.run_summary = run_summary
        self.hgt_summary = hgt_summary
        self.cluster_threshold = cluster_threshold
        self.clonal_clusters = self.get_clonal_clusters(self.cluster_threshold)
        self.full_run_path = config.run_path / f'{self.species_name}__pairwise_runs.pkl'
        with open(self.full_run_path, 'rb') as f:
            self.all_runs = pickle.load(f)

    def get_ani_mat(self, pops=None):
        ani_mat = long_form_to_symmetric(self.hgt_summary, row_name='genome1', col_name='genome2', val_name='ani',
                                         diagonal_fill=1)
        return ani_mat
    
    def get_percid_mat(self):
        id_mat = long_form_to_symmetric(self.hgt_summary, row_name='genome1', col_name='genome2', val_name='perc_id',
                                        diagonal_fill=1)
        return id_mat
    
    def get_syn_div_mat(self):
        div_mat = long_form_to_symmetric(self.run_summary, row_name='genome1', col_name='genome2', val_name='div')
        return div_mat
    
    def get_hgt_mags(self):
        return list(set(self.hgt_summary['genome1'].unique()) | set(self.hgt_summary['genome2'].unique()))
    
    def get_all_mags(self):
        return list(set(self.run_summary['genome1'].unique()) | set(self.run_summary['genome2'].unique()))
    
    def get_clonal_clusters(self, threshold=config.clonal_cluster_pi_threshold):
        """
        Cluster MAGs by percent identical genes
        Requires dRep pairwise results to contain this species
        Returns a DataFrame with columns 'genome' and 'cluster'
        """
        if len(self.hgt_summary) == 0:
            raise ValueError('No HGT results found for this species')

        symmetric_table = long_form_to_symmetric(self.hgt_summary, row_name='genome1', col_name='genome2',
                                                 val_name='perc_id', diagonal_fill=1)
        # convert to perc different so that it's a distance matrix
        symmetric_table = 1 - symmetric_table

        clusters = cluster_close_pairs(symmetric_table, 1-threshold)
        return clusters
    
    def get_clades(self, allowed_pops=None):
        """
        Cluster MAGs by synonymous divergence into two 
        """
        ani_mat = 1 - self.get_syn_div_mat()
        if allowed_pops is not None:
            row_mask = self.metadata.filter_mags_by_pops(ani_mat.index, allowed_pops)
            col_mask = self.metadata.filter_mags_by_pops(ani_mat.columns, allowed_pops)
            ani_mat = ani_mat.loc[row_mask, col_mask]
        if ani_mat.shape[0] < 2:
            return None
        Z = linkage_clustering(1-ani_mat)
        max_clusters = 2
        clusters = fcluster(Z, t=max_clusters, criterion='maxclust')
        mags1 = ani_mat.index[clusters == 1]
        mags2 = ani_mat.index[clusters == 2]
        return mags1, mags2
    
    def clade_statistics(self, mags1, mags2):
        # TODO: port over the codes for clade differentiation statistics
        pass
    
    def set_nonclonal_mags(self, passed_mags=None):
        if passed_mags is None:
            self.nonclonal_mags = sample_genome_per_cluster(self.clonal_clusters)
        else:
            self.nonclonal_mags = passed_mags

    def get_nonclonal_mags(self):
        if not hasattr(self, 'nonclonal_mags'):
            self.set_nonclonal_mags()
        return self.nonclonal_mags
    
    def get_nonclonal_pairwise_summary(self, summary):
        """
        Filter the species pairwise comparison by clonal cluster
        Summary df needs to have columns 'genome1' and 'genome2'
        """
        return filter_df_by_samples(summary, self.get_nonclonal_mags())
    
    def get_pair_full_runs(self, mag1, mag2):
        if (mag1, mag2) in self.all_runs:
            return self.all_runs[(mag1, mag2)]
        elif (mag2, mag1) in self.all_runs:
            return self.all_runs[(mag2, mag1)]
        else:
            raise ValueError(f'Pairwise run not found for {mag1} and {mag2}')
    
    def get_ANI_dist(self, pops):
        if pops is None:
            return self.hgt_summary['ani'].values
        mask = self.hgt_summary['study_x'].isin(pops) & self.hgt_summary['study_y'].isin(pops)
        return self.hgt_summary[mask]['ani'].values
    
    def get_ANI_dist_by_mags(self, mags):
        mask = self.hgt_summary['genome1'].isin(mags) & self.hgt_summary['genome2'].isin(mags)
        return self.hgt_summary[mask]['ani'].values
    
    def get_ANI_dist_between_mags(self, mags1, mags2):
        mask1 = self.hgt_summary['genome1'].isin(mags1) & self.hgt_summary['genome2'].isin(mags2)
        mask2 = self.hgt_summary['genome1'].isin(mags2) & self.hgt_summary['genome2'].isin(mags1)
        return self.hgt_summary[mask1 | mask2]['ani'].values

class PairwiseHelper:
    """
    A class to hold all pairwise comparisons for a databatch; mostly for
    creating species-level helper
    TODO: eventually add the identical gene and genome ANI statistics
    TODO: might need to compute identical block myself because drep cannot finish every pair
    """
    def __init__(self, databatch):
        self.databatch = databatch
        self.metadata = metadata_utils.MetadataHelper(databatch)
        self.pairwise_summary = self.load_pairwise_summary()
        self.hgt_summary = load_hgt_res()

    def get_species_list(self):
        return self.pairwise_summary['species'].unique()

    def load_pairwise_summary(self):
        pairwise_summary_path = config.run_path / f'{self.databatch}_annotated.csv'
        return pd.read_csv(pairwise_summary_path)

    def get_species_pairwise_summary(self, species_name):
        species_res = self.pairwise_summary[self.pairwise_summary['species'] == species_name]
        return species_res
    
    def get_species_helper(self, species_name, cluster_threshold=config.clonal_cluster_pi_threshold):
        species_hgt = self.hgt_summary[self.hgt_summary['species'] == species_name]
        species_run = self.get_species_pairwise_summary(species_name)
        return SpeciesPairwiseHelper(species_name, species_run, species_hgt, metadata=self.metadata, cluster_threshold=cluster_threshold)

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
