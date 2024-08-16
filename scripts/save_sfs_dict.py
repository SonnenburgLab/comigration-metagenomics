import logging
import pandas as pd

from utils import snv_utils, pairwise_utils
import config

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)

all_mag_df = []
for species_name in pairwise_helper.metadata.get_species_list():
    logging.info(f'Processing {species_name}')
    species_helper = pairwise_helper.get_species_helper(species_name)
    passed_mags = species_helper.get_nonclonal_mags()
    mag_df = pd.DataFrame(index=passed_mags, columns=['species'])
    mag_df['species'] = species_name
    all_mag_df.append(mag_df)
    species_snv_helper = snv_utils.SNVHelper(species_name, data_batch=config.databatch)
    species_snv_helper.save_dadi_data_dict(pops=pairwise_helper.metadata.get_all_pops(), allowed_mags=passed_mags)

all_mag_df = pd.concat(all_mag_df)
all_mag_df.to_csv(config.sfs_path / '{}_nonclonal_mags.csv'.format(config.databatch), index=True)