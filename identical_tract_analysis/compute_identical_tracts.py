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

for species in species_list:
    save_path = config.run_path / f'{species}__pairwise_runs.pkl'
    if save_path.exists():
        logging.info(f'Skipping {species}')
        continue
    snv_helper = snv_utils.PairwiseSNVHelper(species_name=species, data_batch=data_batch)
    all_pairs = snv_helper.get_all_genome_pairs()
    logging.info(f'Processing {species}: {len(all_pairs)} pairs')

    all_runs = snv_helper.compute_runs_for_pairs(all_pairs)
    with open(save_path, 'wb') as f:
        pickle.dump(all_runs, f)

    # no longer using identical fraction so this is commented out
    # all_frac_ids = snv_helper.compute_frac_id_for_pairs(all_pairs)
    # save_path = config.identical_fraction_path / f'{species}__identical_fraction.tsv'
    # all_frac_ids.to_csv(save_path, sep='\t', index=False)