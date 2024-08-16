import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
# from matplotlib.backends.backend_pdf import PdfPages

import moments
from utils import moments_utils, metadata_utils, snv_utils
import config


# TODO: set up a yml file for all these parameters
# data related parameters
data_batch = config.databatch
min_sample_size = 30
mask_singletons = True
# focal_pops = ['Hadza', 'Tsimane']
# focal_pops = ['China', 'MetaHIT']
# focal_pops = ['China', 'HMP']
focal_pops = ['MetaHIT', 'HMP']
# focal_pops = ['Nepal', 'MetaHIT']

# set up model related parameters
model = moments.Demographics2D.split_mig
model_name = 'split_mig'
# for split_mig: nu1, nu2, T, m
param_names = ['nu1', 'nu2', 'T', 'm']
p_guess = [2, 2, .1, 1]
lower_bound = [1e-3, 1e-3, 1e-3, 1e-3]
upper_bound = [100, 100, 20, 10]

# model = moments_utils.split_no_mig
# model_name = 'split_no_mig'
# param_names = ['nu1', 'nu2', 'T']
# p_guess = [2, 2, .1]
# lower_bound = [1e-3, 1e-3, 1e-3]
# upper_bound = [100, 100, 20]

# set up output files
# output_pdf = 'figs/moments_split_mig_sfs.pdf'
output_path = Path('moments_out')
output_path.mkdir(exist_ok=True)
output_file = output_path / f'{data_batch}__{model_name}__{focal_pops[0]}__{focal_pops[1]}.csv'

# preparing the dataframe
columns = ['species'] + param_names + ['theta', 'log_likelihood'] \
    + [f'uncert_{param}' for param in param_names] + ['uncert_theta']\
    + ['syn_genome_length', 'num_sites_passing_proj', 'num_syn_snps', 'num_focal_snvs', 'num_snps_after_projection'] \
    + [f'num_{pop}' for pop in focal_pops] \
    + [f'proj_{pop}' for pop in focal_pops]
# result_df = pd.DataFrame(columns=columns)

# write header to file
output_file = open(output_file, 'w')
output_file.write(','.join(columns) + '\n')

metahelper = metadata_utils.MetadataHelper(data_batch=data_batch)
metahelper.species = metahelper.species[~metahelper.species.index.duplicated()]
logging.info(f"Processing total {len(metahelper.species)} species after dropping duplicate rows")
full_species_list = metahelper.get_species_list()

for species in full_species_list:
    logging.info(f"Processing {species}")
    mag_counts = metahelper.get_mag_counts(species)[focal_pops]
    if min(mag_counts) < min_sample_size:
        logging.info(f"Skipping {species} due to low sample size")
        continue

    res = moments_utils.load_SFS(species, metahelper, focal_pops=focal_pops)
    if res is None:
        logging.info(f"Skipping {species} due to insufficient data after projection")
        continue
    data, snv_stats = res
    logging.info(f"Loaded SFS for {species}; {snv_stats}")

    opt_params, theta, uncerts, ll = moments_utils.fit_model(
        data=data, model_func=model, p_guess=p_guess, lower_bound=lower_bound, 
        upper_bound=upper_bound, mask_singletons=mask_singletons)
    
    logging.info(f"Optimal parameters for {species}: {opt_params}")
    logging.info(f"Theta for {species}: {theta}")
    logging.info(f"Uncertainties for {species}: {uncerts}")
    logging.info(f"Log likelihood for {species}: {ll}")

    snv_helper = snv_utils.SNVHelper(species_name=species, data_batch=data_batch)
    syn_genome_length = snv_helper.get_4D_core_genome_length()
    num_syn_snps = snv_helper.get_num_syn_snps()
    assert(snv_stats['num_total_snvs'] == num_syn_snps)
    num_focal_snvs = snv_stats['num_focal_snvs']
    num_snps_after_projection = snv_stats['num_proj_snvs']

    # TODO: added on 2024/05/29. Calculates number of sites that pass the projection threshold
    coverage = snv_helper.get_coverage()
    core_4D_sites = snv_helper.get_4D_core_indices()
    coverage_mat = coverage.loc[core_4D_sites]
    passed_sites = []
    for i, pop in enumerate(focal_pops):
        pop_mask = snv_helper.get_population_mask(pop)
        pop_coverage = coverage_mat.loc[:, pop_mask]
        proj = snv_stats[f'proj_{pop}']
        passed_sites.append(pop_coverage.sum(axis=1) >= proj)
    num_passed = np.sum(passed_sites[0] & passed_sites[1])

    line_vals = [species] + opt_params.tolist() + [theta, ll]  \
                + uncerts.tolist() \
                + [syn_genome_length, num_passed, num_syn_snps, num_focal_snvs, num_snps_after_projection] \
                + mag_counts.values.tolist() + [snv_stats[f'proj_{pop}'] for pop in focal_pops]

    logging.info(f"Writing {species} to file")
    output_file.write(','.join(map(str, line_vals)) + '\n')
    output_file.flush()
output_file.close()
# result_df.to_csv(output_file, index=False)

