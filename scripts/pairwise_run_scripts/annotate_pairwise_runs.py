"""
03/03/2025 -- moved most of the function to prepare_max_run_df.py -- probably can delete this script
"""
import pandas as pd
import pickle
import numpy as np

from utils import metadata_utils, snv_utils
import config

full_res = pd.read_csv(config.run_path / f'{config.databatch}_all_species.csv')

annotated_res = []
for species_name, species_res in full_res.groupby('species'):
    print(species_name)
    save_path = config.run_path / f'{species_name}__pairwise_runs.pkl'
    with open(save_path, 'rb') as f:
        all_runs = pickle.load(f)

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
    
    species_res.set_index(['genome1', 'genome2'], inplace=True)
    species_res.loc[comps, 'num_snvs'] = num_snvs
    species_res.loc[comps, 'core_len'] = core_len
    species_res.reset_index(inplace=True)

    species_res['div'] = species_res['num_snvs'] / species_res['core_len']
    annotated_res.append(species_res)

annotated_res = pd.concat(annotated_res)
annotated_res.to_csv(config.run_path / f'{config.databatch}_annotated.csv', index=False)