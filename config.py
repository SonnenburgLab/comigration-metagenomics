import os
from pathlib import Path
import seaborn as sns

data_base_path = Path('/Volumes/Botein/comigration_data/')
snv_catalog_path = data_base_path / 'snv_catalogs'
# sfs are saved in the "data dict" format specified by dadi, but also usable by moments
sfs_path = data_base_path / 'sfs_data'

# cache snv catalog into a compressed dataframe
cache_snvs_path = data_base_path / 'cache_snvs'
cache_format = 'feather'

dadi_dat_path = data_base_path / 'dadi_data'
vcf_path = data_base_path / 'vcf_data'

# IBS-run analysis path
run_path = data_base_path / 'IBS_runs'

# identical fraction path
identical_fraction_path = data_base_path / 'identical_fraction'

# dRep path
# drep_path = Path('/Users/Device6/Documents/Research/bgoodlab/microbiome_codiv/dRep_analysis')
hgt_res_path = Path('/Users/Device6/Documents/Research/bgoodlab/microbiome_codiv/dRep_analysis/Global_HGT_v7.1.csv')

# close pair thresholds
close_pair_div_ratio = 10
clonal_cluster_pi_threshold = 0.2 # larger than 20% identical gene probably have clonal background
close_pair_pi_threshold = 0.75 # similar to Liu Good 2024

# color pallette for comparing populations
all_comps = ['Tsimane-Tsimane', 'Hadza-Tsimane', 'Hadza-Hadza', 'China-HMP', 'China-China', 'HMP-HMP']
color_dict = dict(zip(all_comps, sns.color_palette('pastel', n_colors=len(all_comps))))

# default databatch
databatch = '240714'