import pandas as pd

import config

class MetadataHelper:
    def __init__(self, data_batch):
        self.data_batch = data_batch
        self.mag_path = config.snv_catalog_path / f'{data_batch}_snv_catalog_mag_metadata.tsv'
        self.species_path = config.snv_catalog_path / f'{data_batch}_snv_catalog_species_metadata.tsv'

        self.mag = pd.read_csv(self.mag_path, sep='\t').set_index('genome')

        species_metadata = pd.read_csv(self.species_path, sep='\t')
        species_metadata['species_names'] = species_metadata['rep_genome'].apply(lambda x: x.split('__')[0])
        species_metadata.set_index('species_names', inplace=True)
        species_metadata.sort_index(inplace=True)
        self.species = species_metadata

        self.mag_to_pop = self.mag['study'].to_dict()

    def get_species_list(self):
        return self.species.index.tolist()

    def get_mag_pop(self, mag_name):
        return self.mag_to_pop.get(mag_name, None)
    
    def get_all_pops(self):
        return self.mag['study'].unique()