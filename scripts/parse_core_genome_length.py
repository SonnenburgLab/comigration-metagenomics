
import numpy as np
import pandas as pd
import os
from utils import snv_utils, dadi_utils
import config

data_batch = '240318_all_species'
species_folder = os.path.join(config.snv_catalog_path, data_batch)

species_list = [species for species in os.listdir(species_folder) if not species.startswith('.')]
species_list = sorted(species_list)
metadata_path = os.path.join(config.snv_catalog_path, '{}__metadata.csv'.format(data_batch))

# generate a metadata along the way
metadata = pd.DataFrame(species_list)
metadata.columns = ['Species']
metadata['Core Genome Len'] = 0
metadata['Core Genome 4D Len'] = 0
metadata.set_index('Species', inplace=True)

for species, _ in metadata.iterrows():
    print(species)

    snv_file = os.path.join(species_folder, species, 'output', '{}.catalog.noAuto.wtRef.tsv'.format(species))
    degeneracy_file = os.path.join(species_folder, species, 'output', '{}_4D_sites.tsv'.format(species))
    core_gene_file = os.path.join(species_folder, species, 'output', '{}_core_genome_mask.tsv'.format(species))

    num_core, num_4D_core = snv_utils.get_core_genome_length(degeneracy_file, core_gene_file)

    metadata.loc[species, 'Core Genome Len'] = num_core
    metadata.loc[species, 'Core Genome 4D Len'] = num_4D_core


metadata.to_csv(metadata_path)