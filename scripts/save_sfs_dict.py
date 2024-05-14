import logging

from utils import metadata_utils, snv_utils

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
# Load the metadata
metadata = metadata_utils.MetadataHelper(data_batch='240509')
pops = metadata.get_all_pops()
for species in metadata.get_species_list():
    logging.info(f'Processing {species}')
    snv_helper = snv_utils.SNVHelper(species_name=species, 
                                     data_batch='240509', cache_snvs=True)

    snv_helper.save_dadi_data_dict(pops=pops)