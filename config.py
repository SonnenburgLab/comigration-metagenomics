from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
project_path = CURRENT_DIR
data_base_path = CURRENT_DIR / 'example'

snv_catalog_path = data_base_path / 'snv_catalogs'
# sfs are saved in the "data dict" format specified by dadi, but also usable by moments
sfs_path = data_base_path / 'sfs_data'

# cache snv catalog into a compressed dataframe
cache_snvs_path = data_base_path / 'cache_snvs'
cache_format = 'feather'

# IBS-run analysis path
run_path = data_base_path / 'IBS_tracts'

# identical fraction path
identical_fraction_path = data_base_path / 'identical_fraction'

# moments ananlysis path
moments_path = project_path / 'moments_analysis'

# ibs tract analysis path
ibs_analysis_path = project_path / 'identical_tract_analysis'

# intermediate data path
intermediate_data_path = project_path / 'temp_dat'

# close pair thresholds
clonal_cluster_pi_threshold = 0.2 # for clustering clonal MAGs
fig2_perc_id_threshold = 0.1

# number of pairs to compute L99
fig3_pair_threshold = 200

# default databatch; the date-id for identifying the SNV catalog and SFS
databatch = 'test'
sfs_batch = f'{databatch}_full'

# cphmm result path
cphmm_result_path = data_base_path / 'cphmm_results'

# mutation clock
mut_rate = 1e-9
gen_per_day = 1
day_per_year = 365.25
mut_per_year = mut_rate * gen_per_day * day_per_year
