import pandas as pd
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utils import metadata_utils, snv_utils
import config

metadata = metadata_utils.MetadataHelper(config.databatch)
syn_genome_lens = pd.Series(index=metadata.species.index)
for species in metadata.species.index:
    print(f'Processing {species}')
    # if first time running this, will create the cache
    snv_helper = snv_utils.SNVHelper(species, data_batch=config.databatch, cache_snvs=True)
    syn_genome_len = snv_helper.get_4D_core_genome_length()
    syn_genome_lens.loc[species] = syn_genome_len

syn_genome_lens.to_csv(config.cache_snvs_path / config.databatch / 'syn_genome_lens.csv')