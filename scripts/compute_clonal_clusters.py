import pandas as pd
import numpy as np
import os

from utils import pairwise_utils
import config

pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)

metadata = pairwise_helper.metadata

all_dfs = []
for species in metadata.get_species_list():
    species_helper = pairwise_helper.get_species_helper(species)
    try:
        cluster_df = species_helper.get_clonal_clusters()
    except ValueError:
        print(f'No clonal clusters for {species}')
        continue
    cluster_df['pop'] = list(map(metadata.get_mag_pop, cluster_df.index))
    cluster_df['species_name'] = species 
    all_dfs.append(cluster_df)

full_clonal_df = pd.concat(all_dfs)

# prepare summary table
cluster_summary_df = pd.DataFrame(index=full_clonal_df['species_name'].unique(), columns=['num_clusters', 'largest_cluster', 'num_mags', 'num_in_cluster'])

for species, clonal_df in full_clonal_df.groupby('species_name'):
    cluster_summary_df.loc[species, 'num_mags'] = clonal_df.shape[0]

    cluster_sizes = clonal_df.groupby('cluster').size()
    cluster_summary_df.loc[species, 'num_unrelated'] = len(cluster_sizes)
    cluster_sizes = cluster_sizes[cluster_sizes > 1]
    cluster_summary_df.loc[species, 'num_clusters'] = len(cluster_sizes)
    if cluster_sizes.shape[0] == 0:
        cluster_summary_df.loc[species, 'largest_cluster'] = 0
    else:
        cluster_summary_df.loc[species, 'largest_cluster'] = cluster_sizes.max()

    cluster_summary_df.loc[species, 'num_in_cluster'] = cluster_sizes.sum()

    fractions = cluster_sizes / cluster_sizes.sum()
    entropy = -np.sum(fractions * np.log(fractions))
    cluster_summary_df.loc[species, 'cluster_entropy'] = entropy

    for pop, pop_df in clonal_df.groupby('pop'):
        num_unrelated = pop_df.groupby('cluster').size().shape[0]
        cluster_summary_df.loc[species, f'num_{pop}_unrelated'] = num_unrelated

cluster_summary_df['frac_in_cluster'] = cluster_summary_df['num_in_cluster'] / cluster_summary_df['num_mags']
cluster_summary_df['frac_largest'] = cluster_summary_df['largest_cluster'] / cluster_summary_df['num_mags']
cluster_summary_df['mean_cluster_frac'] = cluster_summary_df['frac_in_cluster'] / cluster_summary_df['num_clusters']
cluster_summary_df['frac_unrelated'] = cluster_summary_df['num_unrelated'] / cluster_summary_df['num_mags']

cluster_summary_df['HT_unrelated'] = cluster_summary_df['num_Hadza_unrelated'] + cluster_summary_df['num_Tsimane_unrelated']
cluster_summary_df.sort_values('HT_unrelated', ascending=False, inplace=True)

if not os.path.exists('./dat'):
    os.makedirs('./dat')
cluster_summary_df.to_csv('./dat/240816_clonal_cluster_summary.tsv', sep='\t')