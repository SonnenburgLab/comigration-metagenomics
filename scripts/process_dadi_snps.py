"""
TODO: can probably safely delete this script 05/14/2024
"""
import numpy as np
import pandas as pd
import os
from utils import snv_utils, dadi_utils
import config

data_batch = '240326_with_industrial'
dadi_dat_path = os.path.join(config.dadi_dat_path, data_batch)
if not os.path.exists(dadi_dat_path):
    os.mkdir(dadi_dat_path)
species_folder = os.path.join(config.snv_catalog_path, data_batch)
pop_names = ['Hadza', 'Tsimane', 'HMP', 'MetaHIT']

species_list = [species for species in os.listdir(species_folder) if not species.startswith('.')]
species_list = sorted(species_list)
metadata_path = os.path.join(config.dadi_dat_path, '{}__metadata.csv'.format(data_batch))

# generate a metadata along the way
metadata = pd.DataFrame(species_list)
metadata.columns = ['Species']
metadata['Core Genome Len'] = 0
metadata['Core Genome 4D Len'] = 0
metadata['Num Hadza'] = 0
metadata['Num Tsimane'] = 0
metadata['Num SNVs'] = 0
# these two fields will change depending on the later batch; might need to automate
metadata['Num HMP'] = 0
metadata['Num MetaHIT'] = 0
metadata.set_index('Species', inplace=True)

# these are strings in the accession that can be used to identify a population
pop_ids = list(map(snv_utils.get_population_accession_identifier, pop_names))

for species, _ in metadata.iterrows():
    print(species)

    snv_file = os.path.join(species_folder, species, 'output', '{}.catalog.noAuto.wtRef.tsv'.format(species))
    degeneracy_file = os.path.join(species_folder, species, 'output', '{}_4D_sites.tsv'.format(species))
    core_gene_file = os.path.join(species_folder, species, 'output', '{}_core_genome_mask.tsv'.format(species))

    num_core, num_4D_core = snv_utils.get_core_genome_length(degeneracy_file, core_gene_file)
    if num_core == 0:
        print('No core genome for', species)
        continue

    syn_snvs, _  = snv_utils.load_and_filter_snv_catalog(snv_file, degeneracy_file, core_gene_file)

    pop_snvs = snv_utils.split_snvs(syn_snvs, pops=pop_ids)
    # dadi_file = os.path.join(dadi_dat_path, species + '.snps.txt')
    # dadi_utils.write_dadi_input(dadi_file, pop_snvs, pops=['Hadza', 'Tsimane', 'HMP', 'MetaHIT'])

    # write metadata
    metadata.loc[species, 'Core Genome Len'] = num_core
    metadata.loc[species, 'Core Genome 4D Len'] = num_4D_core
    metadata.loc[species, 'Num Hadza'] = pop_snvs[0].shape[1]
    metadata.loc[species, 'Num Tsimane'] = pop_snvs[1].shape[1]
    metadata.loc[species, 'Num HMP'] = pop_snvs[2].shape[1]
    metadata.loc[species, 'Num MetaHIT'] = pop_snvs[3].shape[1]
    metadata.loc[species, 'Num SNVs'] = pop_snvs[0].shape[0]

metadata.to_csv(metadata_path)