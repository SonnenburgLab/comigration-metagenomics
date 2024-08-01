"""
Compute the runs for all pairwise comparisons of genomes in the dataset.

Takes a few hours and O(10) GB of memory.
"""
import logging
import pickle
import pandas as pd
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

from utils import metadata_utils, snv_utils
import config

data_batch = config.databatch
# Load the metadata
metadata = metadata_utils.MetadataHelper(data_batch=data_batch)

species_list = metadata.get_species_list()
# on 24.06.27: trying to compute the new species first
# old_metadata = metadata_utils.MetadataHelper(data_batch='240509')
# old_species_list = old_metadata.get_species_list()
# species_list = [s for s in species_list if s not in old_species_list]

# on 24.06.27: trying to only compute missed species
# full_res = pd.read_csv(config.run_path / '240627_all_runs_all_species_annotated.csv')
# old_species = full_res['species'].unique()
# species_list = [s for s in species_list if s not in old_species]
# logging.info(f'Processing {len(species_list)} species')

for species in species_list:
    snv_helper = snv_utils.PairwiseSNVHelper(species_name=species, data_batch=data_batch)
    all_pairs = snv_helper.get_all_genome_pairs()
    logging.info(f'Processing {species}: {len(all_pairs)} pairs')

    # all_runs = snv_helper.compute_runs_for_pairs(all_pairs)
    # save_path = config.run_path / f'{species}__pairwise_runs.pkl'
    # with open(save_path, 'wb') as f:
    #     pickle.dump(all_runs, f)

    all_frac_ids = snv_helper.compute_frac_id_for_pairs(all_pairs)
    save_path = config.identical_fraction_path / f'{species}__identical_fraction.tsv'
    all_frac_ids.to_csv(save_path, sep='\t', index=False)