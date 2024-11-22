import pandas as pd

from utils import metadata_utils, pairwise_utils
import config

def get_filtered_runs(species_helper, perc_id_threshold=0.1):
    df1 = species_helper.run_summary.copy()
    df2 = species_helper.hgt_summary.copy()
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


metadata = metadata_utils.MetadataHelper(data_batch=config.databatch)

species_list = metadata.get_species_list()

pairwise_helper = pairwise_utils.PairwiseHelper(config.databatch)

percid_threshold = 0.8

summaries = []
for species in species_list:
    print(species)
    species_helper = pairwise_helper.get_species_helper(species, cluster_threshold=percid_threshold)
    if species_helper.run_summary.shape[0] == 0:
        print(f'No run data for {species}')
        continue
    dedup_summary = get_filtered_runs(species_helper, perc_id_threshold=percid_threshold)
    if dedup_summary.shape[0] == 0:
        print(f'No non-clonal pair for {species}')
        continue
    # might need to add include_groups=False depending on pandas version
    num_comps = dedup_summary.groupby('comp', group_keys=False).apply(len)
    score = dedup_summary.groupby('comp', group_keys=False).apply(pairwise_utils.compute_L99)

    # combine the two
    comp_summary = pd.concat([num_comps, score], axis=1)
    comp_summary.columns = ['num_comps', '99_perc_length']
    comp_summary['implied_years'] = pairwise_utils.length_to_years(comp_summary['99_perc_length'])
    comp_summary.reset_index(inplace=True)
    comp_summary['species'] = species
    summaries.append(comp_summary)

summary_df = pd.concat(summaries)
summary_df.to_csv(config.intermediate_data_path / f'241016__run_year_summary__percid={percid_threshold}__.tsv', index=False, sep='\t')