import pandas as pd
import numpy as np
import pickle

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils import metadata_utils
import config

metadata = metadata_utils.MetadataHelper(data_batch=config.databatch)

all_results = []
problem_species = []

def prepare_max_run_df(runs):
    dat = [(pair[0], pair[1], max(runs[pair])) for pair in runs]
    return pd.DataFrame(dat, columns=['genome1', 'genome2', 'max_run'])

def sort_and_join(row):
    sorted_values = sorted([row['pop1'], row['pop2']])
    return '-'.join(sorted_values)

def annotate_df(df):
    df['pop1'] = df['genome1'].apply(lambda x: metadata.get_mag_pop(x))
    df['pop2'] = df['genome2'].apply(lambda x: metadata.get_mag_pop(x))
    df['comp'] = df.apply(sort_and_join, axis=1)
    return df

def compute_div_and_len(all_runs):
    # count number of runs and total length of runs; this can be translated to core length and divergence
    num_runs = []
    total_len = []
    comps = []
    for comp, runs in all_runs.items():
        comps.append(comp)
        num_runs.append(len(runs))
        total_len.append(sum(runs))

    num_snvs = np.array(num_runs) - 1
    core_len = np.array(total_len) + num_snvs
    return num_snvs, core_len, comps

if __name__ == '__main__':
    species_list = metadata.get_species_list()

    for species in species_list:
        save_path = config.run_path / config.databatch / f'{species}__pairwise_runs.pkl'
        with open(save_path, 'rb') as f:
            all_runs = pickle.load(f)
        try:
            max_runs = prepare_max_run_df(all_runs)
        except ValueError as e:
            print(e)
            print(f'Error for {species}')
            problem_species.append(species)
            continue
        max_runs['species'] = species

        max_runs = annotate_df(max_runs)

        # annotate with divergence and core length information
        num_snvs, core_len, comps = compute_div_and_len(all_runs)
        max_runs.set_index(['genome1', 'genome2'], inplace=True)
        max_runs.loc[comps, 'num_snvs'] = num_snvs
        max_runs.loc[comps, 'core_len'] = core_len
        max_runs['div'] = max_runs['num_snvs'] / max_runs['core_len']

        all_results.append(max_runs)

    full_res = pd.concat(all_results).reset_index()
    full_res.to_csv(config.run_path / f'{config.databatch}_annotated_max_runs.csv', index=False)