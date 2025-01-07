from pathlib import Path

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

# files on my laptop
project_path = Path('/Users/Device6/Documents/Research/bgoodlab/microbiome_codiv/comigration_metagenomics')
# dRep path
hgt_res_path = Path('/Users/Device6/Documents/Research/bgoodlab/microbiome_codiv/dRep_analysis/Global_HGT_v7.2.csv')
# figure path
fig_path = project_path / 'figs'
# intermediate data path
intermediate_data_path = project_path / 'dat'

# close pair thresholds
close_pair_div_ratio = 10
clonal_cluster_pi_threshold = 0.2 # larger than 20% identical gene probably have clonal background
close_pair_pi_threshold = 0.75 # similar to Liu Good 2024
fig2_perc_id_threshold = 0.1

# number of pairs to compute L99
fig3_pair_threshold = 200

# color palette for comparing populations
all_comps = ['Tsimane-Tsimane', 'Hadza-Tsimane', 'Hadza-Hadza', 'China-HMP', 'China-China', 'HMP-HMP']
colors = ['#a1c9f4', '#ffb482', '#8de5a1', '#ff9f9b', '#d0bbff', '#debb9b']
color_dict = dict(zip(all_comps, colors))
# color palette for 7 populations
pop_colors = ['#a1c9f4', '#ffb482', '#8de5a1', '#ff9f9b', '#d0bbff', '#debb9b', '#fffea3']

# default databatch
databatch = '240714'

# cphmm result path
cphmm_res_path = Path('/Volumes/Botein/comigration_data/close_pair_inference/')
