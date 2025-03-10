import argparse
import pandas as pd
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Parse command-line argument for databatch
parser = argparse.ArgumentParser()
parser.add_argument('--data_batch', type=str, required=True, help='Databatch identifier')
args = parser.parse_args()
data_batch_arg = args.data_batch

from utils import metadata_utils, snv_utils
import config

metadata = metadata_utils.MetadataHelper(data_batch_arg)
syn_genome_lens = pd.Series(index=metadata.species.index)
for species in metadata.species.index:
    print(f'Processing {species}')
    # if first time running this, will create the cache
    snv_helper = snv_utils.SNVHelper(species, data_batch=data_batch_arg, cache_snvs=True)
    syn_genome_len = snv_helper.get_4D_core_genome_length()
    syn_genome_lens.loc[species] = syn_genome_len

output_path = config.cache_snvs_path / data_batch_arg / 'syn_genome_lens.csv'
syn_genome_lens.to_csv(output_path)