import pandas as pd
import pickle

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

species_list = metadata.get_species_list()

for species in species_list:
    save_path = config.run_path / f'{species}__pairwise_runs.pkl'
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
    all_results.append(max_runs)

full_res = pd.concat(all_results)
full_res = annotate_df(full_res)
full_res.to_csv(config.run_path / f'{config.databatch}_all_species.csv', index=False)