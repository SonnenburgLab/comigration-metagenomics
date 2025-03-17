import pandas as pd
import time

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils import metadata_utils, pairwise_utils
import config

start_time = time.time()

databatch = config.databatch

metadata = metadata_utils.MetadataHelper(data_batch=config.databatch)

species_list = metadata.get_species_list()

pairwise_helper = pairwise_utils.PairwiseHelper(config.databatch)

percid_threshold = config.fig2_perc_id_threshold

if not os.path.exists(config.ibs_analysis_path / 'ibs_dat'):
    os.makedirs(config.ibs_analysis_path / 'ibs_dat')

summaries = []
full_pair_results = []
bootstrap_results = []
for species in species_list:
    print(species)
    species_helper = pairwise_helper.get_species_helper(species, cluster_threshold=percid_threshold)
    if species_helper.run_summary.shape[0] == 0:
        print(f'No run data for {species}')
        continue
    dedup_summary = species_helper.get_filtered_runs(perc_id_threshold=percid_threshold)
    if dedup_summary.shape[0] == 0:
        print(f'No non-clonal pair for {species}')
        continue
    full_pair_results.append(dedup_summary)
    # might need to add include_groups=False depending on pandas version
    num_comps = dedup_summary.groupby('comp', group_keys=False).apply(len)
    score = dedup_summary.groupby('comp', group_keys=False).apply(pairwise_utils.compute_L99)

    bootstrap_df = pairwise_utils.compute_L99_bootstrap(dedup_summary, n_bootstrap=100)
    bootstrap_df['species'] = species
    bootstrap_results.append(bootstrap_df)

    # combine the two
    comp_summary = pd.concat([num_comps, score], axis=1)
    comp_summary.columns = ['num_comps', '99_perc_length']
    comp_summary['l99_in_years'] = pairwise_utils.length_to_years(comp_summary['99_perc_length'])
    comp_summary.reset_index(inplace=True)
    comp_summary['species'] = species
    summaries.append(comp_summary)

pairwise_runs = pd.concat(full_pair_results)
pairwise_runs.to_csv(config.run_path / f'{databatch}__pairwise_max_runs__percid={percid_threshold}__.tsv', index=False, sep='\t')

summary_df = pd.concat(summaries)
summary_df.to_csv(config.ibs_analysis_path / 'ibs_dat' / f'{databatch}__run_year_summary__percid={percid_threshold}__.tsv', index=False, sep='\t')

bootstrap_df = pd.concat(bootstrap_results)
bootstrap_df.to_csv(config.ibs_analysis_path / 'ibs_dat' / f'{databatch}__bootstrap_results__percid={percid_threshold}__.tsv', index=False, sep='\t')

print(f'Finished in {time.time() - start_time:.2f} seconds')