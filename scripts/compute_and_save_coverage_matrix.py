import logging
from utils import snv_utils, metadata_utils
import config

# find species to process from metadata
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
# Load the metadata
data_batch = '240509'
metadata = metadata_utils.MetadataHelper(data_batch=data_batch)

for species in metadata.get_species_list():
    coords_files = config.snv_catalog_path / data_batch / species / 'output' / 'coords_files'
    file_names = list(coords_files.glob('*.coords'))
    logging.info(f'Processing {species}')

    if config.cache_format=='feather':
        cache_filename = config.cache_snvs_path / data_batch / 'coverage' / f'{species}_coverage.feather'
    elif config.cache_format=='parquet':
        cache_filename = config.cache_snvs_path / data_batch / 'coverage' / f'{species}_coverage.parquet'

    if cache_filename.exists():
        logging.info(f'Coverage matrix already exists for {species}. Skipping.')
        continue

    snv_helper = snv_utils.SNVHelper(species_name=species, data_batch=data_batch, cache_snvs=True)
    # save time by not loading the snv matrix
    core_gene_idx = snv_utils.load_core_gene_mask(snv_helper.core_gene_file)
    if core_gene_idx is None:
        logging.info(f'No core gene file found for {species}. Skipping.')
        continue
    res = snv_utils.compute_coverage_matrix_from_coords(file_names, core_gene_idx, snv_helper.genome_names)
    res = res.reset_index()

    if config.cache_format=='feather':
        res.to_feather(config.cache_snvs_path / data_batch / 'coverage' / f'{species}_coverage.feather')
    elif config.cache_format=='parquet':
        res.to_parquet(config.cache_snvs_path / data_batch / 'coverage' / f'{species}_coverage.parquet')