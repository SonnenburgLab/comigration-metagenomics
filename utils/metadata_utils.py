import pandas as pd
import seaborn as sns

import config

class MetadataHelper:
    def __init__(self, data_batch):
        self.data_batch = data_batch
        self.mag_path = config.snv_catalog_path / f'{data_batch}_snv_catalog_mag_metadata.tsv'
        self.species_path = config.snv_catalog_path / f'{data_batch}_snv_catalog_species_metadata.tsv'

        self.mag = pd.read_csv(self.mag_path, sep='\t').set_index('genome')

        species_metadata = pd.read_csv(self.species_path, sep='\t')
        species_metadata['species_names'] = species_metadata['rep_genome'].apply(lambda x: x.split('__')[0])
        # a couple rows have the same species name but different abundances
        species_metadata.drop_duplicates(subset=['species_names'], inplace=True)
        species_metadata.set_index('species_names', inplace=True)
        species_metadata.sort_index(inplace=True)
        self.species = species_metadata

        self.mag_to_pop = self.mag['study'].to_dict()
        self.all_pops = self.mag['study'].unique()
        for pop in self.all_pops:
            # parse str to int
            self.species[pop] = self.species[pop].astype(int)

    def drop_duplicate_species(self):
        self.species = self.species[~self.species.index.duplicated()]

    def get_species_list(self):
        return self.species.index.tolist()
    
    def get_new_species_list(self, old_data_batch):
        old_metadata = MetadataHelper(old_data_batch)
        old_species_list = old_metadata.get_species_list()
        return [s for s in self.get_species_list() if s not in old_species_list]
    
    def get_mag_counts(self, species_name):
        if species_name in self.species.index:
            return self.species.loc[species_name, self.all_pops]
        else:
            raise ValueError(f'{species_name} not found in species metadata')

    def get_mag_pop(self, mag_name):
        return self.mag_to_pop.get(mag_name, None)
    
    def get_all_pops(self):
        return self.all_pops
    
    def get_species_rep_genome(self, species_name):
        return self.species.loc[species_name, 'rep_genome']
    
    def get_mags_pop_count(self, mags):
        pops = [self.get_mag_pop(mag) for mag in mags]
        return pd.Series(pops).value_counts()
    
    def filter_mags_by_pops(self, mags, pops):
        return [mag for mag in mags if self.get_mag_pop(mag) in pops]
    
    def mag_to_pop_colors(self, mags):
        pop_assignments = [self.get_mag_pop(x) for x in mags]
        all_pops = self.get_all_pops()
        palette = sns.color_palette("pastel", len(all_pops))  # You can choose different palettes
        color_dict = {pop: color for pop, color in zip(all_pops, palette)}
        # make paleofeces stand out more
        color_dict['Paleofeces'] = 'black'
        row_colors = [color_dict[pop] for pop in pop_assignments]
        return row_colors