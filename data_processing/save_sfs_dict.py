import logging
import pandas as pd
import argparse

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils import snv_utils, pairwise_utils, metadata_utils
import config

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Save SFS dict for SNV analysis.')
parser.add_argument('--data_batch', required=True, help='Data batch to process')
parser.add_argument('--saving_full', type=lambda s: s.lower() in ['true', '1', 'yes'], default=True,
                    help='Whether to save full version (True) or not (False)')
args = parser.parse_args()

data_batch = args.data_batch
saving_full = args.saving_full

metadata = metadata_utils.MetadataHelper(data_batch=data_batch)
all_mag_df = []
if saving_full:
    output_path = config.sfs_path / f'{data_batch}_full'
    for species_name in metadata.get_species_list():
        logging.info(f'Processing {species_name}')
        species_snv_helper = snv_utils.SNVHelper(species_name, data_batch=data_batch)
        species_snv_helper.save_dadi_data_dict(pops=metadata.get_all_pops(), allowed_mags=None, output_path=output_path)
else:
    pairwise_helper = pairwise_utils.PairwiseHelper(databatch=data_batch)
    output_path = config.sfs_path / data_batch
    for species_name in metadata.get_species_list():
        logging.info(f'Processing {species_name}')
        species_snv_helper = snv_utils.SNVHelper(species_name, data_batch=data_batch)
        species_helper = pairwise_helper.get_species_helper(species_name)
        passed_mags = species_helper.get_nonclonal_mags()
        mag_df = pd.DataFrame(index=passed_mags, columns=['species'])
        mag_df['species'] = species_name
        all_mag_df.append(mag_df)
        species_snv_helper.save_dadi_data_dict(pops=pairwise_helper.metadata.get_all_pops(), allowed_mags=passed_mags, output_path=output_path)
    all_mag_df = pd.concat(all_mag_df)
    all_mag_df.to_csv(config.sfs_path / f'{data_batch}_nonclonal_mags.csv', index=True)