import numpy as np
import pandas as pd
import scipy
import pickle
import os
from scipy.cluster.hierarchy import linkage, fcluster
from utils import metadata_utils

import config

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
def compute_L99(df, val_name='max_run', quantile=0.99):
    l99 = df[val_name].quantile(quantile, interpolation='higher')
    return l99

def compute_L99_resampled(df, val_name='max_run', n_bootstrap=100, n_resample=None):
    if n_resample is None:
        n_resample = len(df)
    bootstrap_scores = []
    for _ in range(n_bootstrap):
        # Resample rows with replacement
        resampled_df = df.sample(n=n_resample, replace=True)
        bootstrap_scores.append(compute_L99(resampled_df, val_name=val_name))
    return bootstrap_scores

def _compute_L99_bootstrap(input_df, n_bootstrap=100):
    bootstrap_scores = []
    for _ in range(n_bootstrap):
        # Resample rows with replacement
        resampled_df = input_df.sample(n=len(input_df), replace=True)
        bootstrap_scores.append(compute_L99(resampled_df))
    return bootstrap_scores

def compute_L99_bootstrap(summary, n_bootstrap=100):
    bootstrap_df = pd.DataFrame()
    for name, grouped in summary.groupby('comp', group_keys=False):
        bootstrap_df[name] = _compute_L99_bootstrap(grouped, n_bootstrap)
    bootstrap_df = bootstrap_df.T.reset_index(names='comp')
    return bootstrap_df


def HGT_score(df, quantile=0.99, ref_length=700):
    """
    Idea: take the distribution of max_runs, find the tail length at a certain quantile,
    then compare that with a reference length to get a score

    For example, 700bp corresponds to a length scale of 1 SNV in 10,000 years
    """
    hgt_len = df['max_run'].quantile(quantile, interpolation='higher')
    hgt_score = np.log10(hgt_len / ref_length)
    return hgt_score

def length_to_years(run, mu=config.mut_rate, gen_per_day=config.gen_per_day):
    mu_per_year = mu * gen_per_day * config.day_per_year
    return 1 / run / mu_per_year

def year_to_length(year, mu=config.mut_rate, gen_per_day=config.gen_per_day):
    mu_per_year = mu * gen_per_day * config.day_per_year
    return 1 / year / mu_per_year

def long_run_length(df, quantile=0.99):
    return df['max_run'].quantile(quantile, interpolation='higher')


"""
Below are some classes for organizing pairwise comparisons
"""
class SpeciesPairwiseHelper:
    def __init__(self, species_name, run_summary, drep_summary, metadata, cluster_threshold=config.clonal_cluster_pi_threshold):
        self.metadata = metadata
        self.species_name = species_name
        self.run_summary = run_summary
        self.drep_summary = drep_summary
        self.cluster_threshold = cluster_threshold
        self.clonal_clusters = self.get_clonal_clusters(self.cluster_threshold)
        self.data_batch = metadata.data_batch
        self.full_run_path = config.run_path / self.data_batch / f'{self.species_name}__pairwise_runs.pkl'
        with open(self.full_run_path, 'rb') as f:
            self.all_runs = pickle.load(f)

    def get_ani_mat(self, pops=None):
        ani_mat = long_form_to_symmetric(self.drep_summary, row_name='genome1', col_name='genome2', val_name='ani',
                                         diagonal_fill=1)
        return ani_mat
    
    def get_percid_mat(self):
        id_mat = long_form_to_symmetric(self.drep_summary, row_name='genome1', col_name='genome2', val_name='perc_id',
                                        diagonal_fill=1)
        return id_mat
    
    def get_syn_div_mat(self):
        div_mat = long_form_to_symmetric(self.run_summary, row_name='genome1', col_name='genome2', val_name='div')
        return div_mat
    
    def get_drep_mags(self):
        return list(set(self.drep_summary['genome1'].unique()) | set(self.drep_summary['genome2'].unique()))
    
    def get_all_mags(self):
        return list(set(self.run_summary['genome1'].unique()) | set(self.run_summary['genome2'].unique()))
    
    def get_clonal_clusters(self, threshold=config.clonal_cluster_pi_threshold):
        """
        Cluster MAGs by percent identical genes
        Requires dRep pairwise results to contain this species
        Returns a DataFrame with columns 'genome' and 'cluster'
        """
        if len(self.drep_summary) == 0:
            raise ValueError('No HGT results found for this species')

        symmetric_table = long_form_to_symmetric(self.drep_summary, row_name='genome1', col_name='genome2',
                                                 val_name='perc_id', diagonal_fill=1)
        # convert to perc different so that it's a distance matrix
        symmetric_table = 1 - symmetric_table

        clusters = cluster_close_pairs(symmetric_table, 1-threshold)
        return clusters
    
    def get_filtered_runs(self, perc_id_threshold=0.1):
        df1 = self.run_summary.copy()
        df2 = self.drep_summary.copy()
        df1['genome_pair'] = df1.apply(lambda row: tuple(sorted([row['genome1'], row['genome2']])), axis=1)
        df2['genome_pair'] = df2.apply(lambda row: tuple(sorted([row['genome1'], row['genome2']])), axis=1)

        # Drop the original genome1 and genome2 (optional)
        df1 = df1.drop(columns=['genome1', 'genome2'])
        df2 = df2.drop(columns=['genome1', 'genome2'])

        # Now join the two dataframes on the 'genome_pair' column
        merged_df = pd.merge(df1, df2, on='genome_pair', how='right')

        merged_df[['genome1', 'genome2']] = pd.DataFrame(merged_df['genome_pair'].tolist(), index=merged_df.index)

        # Drop the 'genome_pair' column (optional)
        merged_df = merged_df.drop(columns=['genome_pair'])

        # Set the new 'genome1' and 'genome2' columns as a MultiIndex
        merged_df.set_index(['genome1', 'genome2'], inplace=True)
        merged_df = merged_df[merged_df['perc_id'] < perc_id_threshold]
        return merged_df

    def filter_ani_by_pops(self, ani_type=None, allowed_pops=None):
        if ani_type=='ANI':
            ani_mat = self.get_ani_mat()
        elif ani_type=='syn_div':
            ani_mat = 1 - self.get_syn_div_mat()
        else:
            # default to synoymous divergence of the core genome
            ani_mat = 1 - self.get_syn_div_mat()
        
        if allowed_pops is not None:
            row_mask = self.metadata.filter_mags_by_pops(ani_mat.index, allowed_pops)
            col_mask = self.metadata.filter_mags_by_pops(ani_mat.columns, allowed_pops)
            ani_mat = ani_mat.loc[row_mask, col_mask]
        return ani_mat
    
    def get_clades(self, ani_type=None, allowed_pops=None):
        """
        Cluster MAGs by synonymous divergence into two 
        """
        ani_mat = self.filter_ani_by_pops(ani_type=ani_type, allowed_pops=allowed_pops)
        if ani_mat.shape[0] < 2:
            return None
        Z = linkage_clustering(1-ani_mat)
        mag_names = ani_mat.index.values

        max_clusters = 2
        distance_column = Z[:, 2]
        # find two clusters with the biggest distance
        clusters = fcluster(Z, t=max_clusters, criterion='maxclust_monocrit', monocrit=distance_column)

        # maxclust only actually produces the same results?
        # clusters = fcluster(Z, t=max_clusters, criterion='maxclust')
        mags1 = mag_names[clusters == 1]
        mags2 = mag_names[clusters == 2]
        return mags1, mags2
    
    def visualize_clades(self, ani_type=None, allowed_pops=None):
        ani_mat = self.filter_ani_by_pops(ani_type=ani_type, allowed_pops=allowed_pops)
        Z = linkage_clustering(1-ani_mat)
        max_clusters = 2
        clusters = fcluster(Z, t=max_clusters, criterion='maxclust')

        colors = np.array(['red', 'blue'])
        col_colors = colors[clusters-1]
        row_colors = self.metadata.mag_to_pop_colors(ani_mat.index)

        import seaborn as sns
        g = sns.clustermap(ani_mat, row_linkage=Z, col_linkage=Z, row_colors=row_colors, col_colors=col_colors, 
                           figsize=(10, 10), cbar_kws=dict(orientation='horizontal'))
        x0, y0, w, h = g.cbar_pos
        g.ax_cbar.set_position([x0, 0.95, g.ax_row_dendrogram.get_position().width * 0.9, 0.02])
        g.ax_cbar.set_title(f'ANI')
        return g

    
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
    
    def get_ANI_dist_by_mags(self, mags):
        mask = self.drep_summary['genome1'].isin(mags) & self.drep_summary['genome2'].isin(mags)
        return self.drep_summary[mask]['ani'].values
    
    def get_ANI_dist_between_mags(self, mags1, mags2):
        mask1 = self.drep_summary['genome1'].isin(mags1) & self.drep_summary['genome2'].isin(mags2)
        mask2 = self.drep_summary['genome1'].isin(mags2) & self.drep_summary['genome2'].isin(mags1)
        return self.drep_summary[mask1 | mask2]['ani'].values
    
    def get_between_pop_L99(self, pop1, pop2, perc_id_threshold=0.1):
        if perc_id_threshold is None:
            summary = self.run_summary
        else:
            summary = self.get_filtered_runs(perc_id_threshold)
        mask1 = (summary['pop1']==pop1) & (summary['pop2']==pop2)
        mask2 = (summary['pop1']==pop2) & (summary['pop2']==pop1)
        filtered = summary[(mask1 | mask2)]
        l99 = compute_L99(filtered)
        num_pairs = len(filtered)
        return l99, num_pairs


class PairwiseHelper:
    """
    A class to hold all pairwise comparisons for a databatch; mostly for
    creating species-level helper
    """
    def __init__(self, databatch):
        self.databatch = databatch
        self.metadata = metadata_utils.MetadataHelper(databatch)
        self.pairwise_summary = self.load_pairwise_summary(databatch)
        self.drep_summary = self.load_drep_res(databatch)

    def get_species_list(self):
        return self.pairwise_summary['species'].unique()

    @staticmethod
    def load_pairwise_summary(databatch):
        filepath = config.run_path / f'{databatch}_annotated_max_runs.csv'
        return pd.read_csv(filepath)

    @staticmethod
    def load_drep_res(databatch):
        drep_res_path = config.drep_res_path / f'{databatch}_drep.csv'
        # loading drep results (Matt O's analysis)
        if not os.path.exists(drep_res_path):
            raise FileNotFoundError(f'No dRep results found at {drep_res_path}')
        drep_res = pd.read_csv(drep_res_path)
        # drop the suffixes
        drep_res['species'] = drep_res['name'].str.replace(r'(_v7|_v8)', '', regex=True)
        drep_res['genome1'] = drep_res['querry'].str.replace('.fa', '')
        drep_res['genome2'] = drep_res['reference'].str.replace('.fa', '')
        return drep_res

    def get_species_pairwise_summary(self, species_name):
        species_res = self.pairwise_summary[self.pairwise_summary['species'] == species_name]
        return species_res
    
    def get_species_helper(self, species_name, cluster_threshold=config.clonal_cluster_pi_threshold):
        species_hgt = self.drep_summary[self.drep_summary['species'] == species_name]
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
        species_drep_res = self.drep_summary[self.drep_summary['species'] == species_name]
        symmetric_table = long_form_to_symmetric(species_drep_res, row_name='genome1', col_name='genome2', val_name='perc_id')
        np.fill_diagonal(symmetric_table.values, 1)
        symmetric_table = 1 - symmetric_table

        passed_mags = cluster_close_pairs(symmetric_table, 1-threshold)
        return passed_mags

    def filter_population(self, species_res, pop):
        return species_res[~((species_res['pop1'] == pop) | (species_res['pop2'] == pop))]
