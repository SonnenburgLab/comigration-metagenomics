import logging
from utils import snv_utils, metadata_utils
import config
import pandas as pd

# find species to process from metadata
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
# Load the metadata
# data_batch = '240626'
data_batch = config.databatch
metadata = metadata_utils.MetadataHelper(data_batch=data_batch)

species_list = metadata.get_species_list()

# on 24.06.27: trying to compute the new species first
# old_metadata = metadata_utils.MetadataHelper(data_batch='240509')
# old_species_list = old_metadata.get_species_list()
# species_list = [s for s in species_list if s not in old_species_list]
# full_res = pd.read_csv(config.run_path / '240627_all_runs_all_species_annotated.csv')
# old_species = full_res['species'].unique()
# species_list = [s for s in species_list if s not in old_species]
# logging.info(f'Processing {len(species_list)} species')

savepath = config.cache_snvs_path / data_batch / 'coverage'
savepath.mkdir(parents=True, exist_ok=True)

for species in species_list:
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
        res.to_feather(savepath / f'{species}_coverage.feather')
    elif config.cache_format=='parquet':
        res.to_parquet(savepath / f'{species}_coverage.parquet')