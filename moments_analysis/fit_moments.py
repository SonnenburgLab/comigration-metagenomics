import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
# from matplotlib.backends.backend_pdf import PdfPages
import moments
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils import moments_utils, metadata_utils, snv_utils
import config

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Fit moments model to SFS data.')
    parser.add_argument('--model_to_fit', type=str, choices=['split_mig', 'split_no_mig'], required=True, help='Model to fit')
    parser.add_argument('--sfs_batch', type=str, required=True, help='SFS batch to use')
    parser.add_argument('--focal_pops', type=str, nargs=2, required=True, help='Pairs of populations to use')
    parser.add_argument('--min_sample_size', type=int, default=30, help='Minimum sample size for each population')
    parser.add_argument('--num_reps', type=int, default=5, help='Number of repetitions')
    parser.add_argument('--mask_singletons', type=bool, default=True, help='Mask singletons')
    args = parser.parse_args()

    #################################################
    # Choose these parameters
    #################################################
    model_to_fit = args.model_to_fit
    sfs_batch = args.sfs_batch
    focal_pops = args.focal_pops
    min_sample_size = args.min_sample_size
    num_reps = args.num_reps
    mask_singletons = args.mask_singletons

    data_batch = sfs_batch.split('_')[0]
    #################################################
    # Set up default parameters for model and data
    #################################################
    # data related parameters
    output_path = config.moments_path / 'moments_dat'
    # sfs_batch = config.sfs_batch
    sfs_dir = config.sfs_path / sfs_batch
    output_path.mkdir(exist_ok=True)

    # set up model related parameters
    if model_to_fit == 'split_mig':
        model = moments.Demographics2D.split_mig
        model_name = 'split_mig'
        # for split_mig: nu1, nu2, T, m
        param_names = ['nu1', 'nu2', 'T', 'm']
        p_guess = [2, 2, .1, 1]
        lower_bound = [1e-3, 1e-3, 1e-3, 1e-3]
        upper_bound = [100, 100, 20, 10]
    elif model_to_fit == 'split_no_mig':
        # same as above but with m fixed to 0
        model = moments_utils.split_no_mig
        model_name = 'split_no_mig'
        param_names = ['nu1', 'nu2', 'T']
        p_guess = [2, 2, .1]
        lower_bound = [1e-3, 1e-3, 1e-3]
        upper_bound = [100, 100, 20]

    # preparing the dataframe for holding results
    columns = ['species', 'rep'] + param_names + ['theta', 'log_likelihood'] \
        + [f'uncert_{param}' for param in param_names] + ['uncert_theta']\
        + ['syn_genome_length', 'num_sites_passing_proj', 'num_syn_snps', 'num_focal_snvs', 'num_snps_after_projection'] \
        + [f'num_{pop}' for pop in focal_pops] \
        + [f'proj_{pop}' for pop in focal_pops]

    # write header to file
    output_file = output_path / f'{sfs_batch}__{model_name}__{focal_pops[0]}__{focal_pops[1]}.csv'
    output_file = open(output_file, 'w')
    output_file.write(','.join(columns) + '\n')

    #################################################
    # Load data and fit model
    #################################################
    metahelper = metadata_utils.MetadataHelper(data_batch=data_batch)
    metahelper.species = metahelper.species[~metahelper.species.index.duplicated()]
    logging.info(f"Processing total {len(metahelper.species)} species after dropping duplicate rows")
    full_species_list = metahelper.get_species_list()

    for species in full_species_list:
        logging.info(f"Processing {species}")
        mag_counts = metahelper.get_mag_counts(species)[focal_pops]
        logging.info(f"Counts for {species}: {mag_counts}")
        if min(mag_counts) < min_sample_size:
            logging.info(f"Skipping {species} due to low sample size")
            continue

        res = moments_utils.load_SFS_projection(species, focal_pops=focal_pops, sfs_folder=sfs_dir)
        if res is None:
            logging.info(f"Skipping {species} due to insufficient data after projection")
            continue
        data, snv_stats = res
        logging.info(f"Loaded SFS for {species}; {snv_stats}")

        # get the syn genome length and number of syn snps
        # important for scaling the inferred time 
        snv_helper = snv_utils.SNVHelper(species_name=species, data_batch=data_batch)
        syn_genome_length = snv_helper.get_4D_core_genome_length()
        num_syn_snps = snv_helper.get_num_syn_snps()
        assert(snv_stats['num_total_snvs'] == num_syn_snps)
        num_focal_snvs = snv_stats['num_focal_snvs']
        num_snps_after_projection = snv_stats['num_proj_snvs']

        # computing the number of sites passing projection requirement
        # basically the sites that are covered by at least the projected number of samples
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

        for rep in range(num_reps):
            logging.info(f"Fitting {species} for rep {rep+1}")
            opt_params, theta, uncerts, ll = moments_utils.fit_model(
                data=data, model_func=model, p_guess=p_guess, lower_bound=lower_bound, 
                upper_bound=upper_bound, mask_singletons=mask_singletons)
            
            logging.info(f"Optimal parameters for {species}: {opt_params}")
            logging.info(f"Theta for {species}: {theta}")
            logging.info(f"Uncertainties for {species}: {uncerts}")
            logging.info(f"Log likelihood for {species}: {ll}")

            line_vals = [species, rep] + opt_params.tolist() + [theta, ll]  \
                        + uncerts.tolist() \
                        + [syn_genome_length, num_passed, num_syn_snps, num_focal_snvs, num_snps_after_projection] \
                        + mag_counts.values.tolist() + [snv_stats[f'proj_{pop}'] for pop in focal_pops]

            logging.info(f"Writing {species} to file")
            output_file.write(','.join(map(str, line_vals)) + '\n')
        output_file.flush()
    output_file.close()

if __name__ == "__main__":
    main()
