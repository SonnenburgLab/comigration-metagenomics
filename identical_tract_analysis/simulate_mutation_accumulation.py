import pandas as pd
import numpy as np
import time

from identical_tract_analysis import identical_run_sim
from utils import metadata_utils, snv_utils, pairwise_utils
import config

def get_passed_species(focal_pops=['Hadza', 'Tsimane']):
    percid_threshold = config.fig2_perc_id_threshold
    run_summary_path = config.intermediate_data_path / f'241016__run_year_summary__percid={percid_threshold}__.tsv'
    run_summary = pd.read_csv(run_summary_path, sep='\t')
    run_summary.set_index(['species', 'comp'], inplace=True)

    run_summary = run_summary[run_summary['num_comps']>config.fig3_pair_threshold]

    passed_species = []
    for species in run_summary.index.get_level_values('species').unique():
        # check if the species has the focal populations
        if_pop1 = (species, '{}-{}'.format(focal_pops[0], focal_pops[0])) in run_summary.index
        if_pop2 = (species, '{}-{}'.format(focal_pops[1], focal_pops[1])) in run_summary.index
        if_pair = (species, '{}-{}'.format(focal_pops[0], focal_pops[1])) in run_summary.index
        if if_pop1 and if_pop2 and if_pair:
            passed_species.append(species)
    return passed_species


def get_passed_species_full_within():
    # let me get some statistics
    percid_threshold = config.fig2_perc_id_threshold
    run_summary_path = config.intermediate_data_path / f'241016__run_year_summary__percid={percid_threshold}__.tsv'
    run_summary_full = pd.read_csv(run_summary_path, sep='\t')
    run_summary_full.set_index(['species', 'comp'], inplace=True)

    run_summary = run_summary_full[run_summary_full['num_comps']>config.fig3_pair_threshold]

    run_summary.reset_index(inplace=True)
    run_summary_full.reset_index(inplace=True)

    all_focal_pops = [['Hadza', 'Tsimane'], ['China', 'HMP'], ['HMP', 'MetaHIT']]
    between_comps = ['-'.join(x) for x in all_focal_pops]
    between_summary = run_summary[run_summary['comp'].isin(between_comps)]

    within_pops = ['{}-{}'.format(pop, pop) for pop in ['Hadza', 'Tsimane', 'China', 'HMP', 'MetaHIT']]
    within_comps = run_summary_full[run_summary_full['comp'].isin(within_pops)]
    within_pairs = within_comps.groupby('species')['num_comps'].sum()

    # find species both in between df and within df
    between_species = between_summary['species'].unique()
    within_species = within_pairs.index.unique()

    common_species = np.intersect1d(between_species, within_species)

    passed_within_pairs = within_pairs.loc[common_species]

    passed_within_pairs = passed_within_pairs[passed_within_pairs.loc[common_species] > 200]
    print("Total {} species and {} pairs".format(len(passed_within_pairs), passed_within_pairs.sum()))
    return passed_within_pairs.index

def rate_to_years(rate, mut_rate=4.08e-10, gen_per_day=1):
    mut_per_year = mut_rate * 365 * gen_per_day
    return rate / mut_per_year

def years_to_rate(years, mut_rate=4.08e-10, gen_per_day=1):
    mut_per_year = mut_rate * 365 * gen_per_day
    return mut_per_year * years

def simulate_one_species(species_helper, focal_pops, years_to_sim, num_pairs=500):
    run_summary = species_helper.get_filtered_runs(perc_id_threshold=config.fig2_perc_id_threshold)
    within_pops = ['{}-{}'.format(pop, pop) for pop in focal_pops]

    within_sum = run_summary[run_summary['comp'].isin(within_pops)]
    # within_runs = within_sum['max_run'].values

    # take some pairs from within populations to simulate mutation accumulation
    if num_pairs < len(within_sum):
        pairs_to_sim = within_sum.sample(num_pairs, replace=False)
    else:
        pairs_to_sim = within_sum

    print("Simulating {} pairs".format(len(pairs_to_sim)))

    pairs_to_sim.reset_index(inplace=True)
    # prepare a dataframe to store the results
    # columns are the simulated max run for each pair
    res_df = pd.DataFrame(columns=['genome1', 'genome2'])
    res_df['genome1'] = pairs_to_sim['genome1']
    res_df['genome2'] = pairs_to_sim['genome2']
    res_df.set_index(['genome1', 'genome2'], inplace=True)

    rates = [years_to_rate(x) for x in years_to_sim]
    for _, row in pairs_to_sim.iterrows():
        mag1 = row['genome1']
        mag2 = row['genome2']
        full_runs = species_helper.get_pair_full_runs(mag1, mag2)
        for i, rate in enumerate(rates):
            sim_runs = identical_run_sim.simulate_mut_accumulation_one_genome(full_runs, rate)
            sim_max = np.max(sim_runs)
            res_df.loc[(mag1, mag2), f'sim_max_run_{i}'] = sim_max
    return res_df


if __name__ == '__main__':
    start_time = time.time()
    pairwise_helper = pairwise_utils.PairwiseHelper(config.databatch)

    # set up simulation parameters
    years_to_sim = np.logspace(3, 5.5, 20)
    num_pairs = 1e6 # no longer subsampling pairs so choose a large number

    passed_within_species = get_passed_species_full_within()
    for species in passed_within_species:
        species_helper = pairwise_helper.get_species_helper(species)
        save_path = config.intermediate_data_path / 'mutation_accumulation' / f'{species}_full_mutation_accumulation.tsv'
        if save_path.exists():
            print(f'{species} already exists. Skip.')
            continue
        res_df = simulate_one_species(species_helper, ['Hadza', 'Tsimane', 'China', 'HMP', 'MetaHIT'], years_to_sim, num_pairs=num_pairs)
        res_df.to_csv(save_path, sep='\t')
        print(f'Simulation for {species} done.')
    print(f'Total time: {time.time() - start_time} seconds')