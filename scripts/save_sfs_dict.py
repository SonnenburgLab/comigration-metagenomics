import logging

from utils import metadata_utils, snv_utils
import config

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

data_batch = config.databatch
# Load the metadata
metadata = metadata_utils.MetadataHelper(data_batch=data_batch)
pops = metadata.get_all_pops()
for species in metadata.get_species_list():
    logging.info(f'Processing {species}')
    snv_helper = snv_utils.SNVHelper(species_name=species, 
                                     data_batch=data_batch, cache_snvs=True)

    snv_helper.save_dadi_data_dict(pops=pops)