import pandas as pd
import numpy as np
import argparse
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import config
from utils import moments_utils

def rescale_time(moment_results):
    # use mutation rate to translate scaled time to years
    moment_results['Tsplit'] = moments_utils.rescale_time(moment_results['T'], moment_results['theta'], moment_results['num_sites_passing_proj'])
    moment_results['Tsplit_uncert'] = moments_utils.rescale_time(moment_results['uncert_T'], moment_results['theta'], moment_results['num_sites_passing_proj'])
    return moment_results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Postprocess moments results.')
    parser.add_argument('--pops', type=str, nargs=2, required=True, help='Pairs of populations to use')
    parser.add_argument('--model_name', type=str, choices=['split_mig', 'split_no_mig'], required=True, help='Model name')
    parser.add_argument('--sfs_batch', type=str, required=True, help='SFS batch to use')
    args = parser.parse_args()

    pops = args.pops
    model_name = args.model_name
    sfs_batch = args.sfs_batch

    inferred_path = config.moments_path / 'moments_dat' / \
        f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}.csv'
    residual_path = config.moments_path / 'moments_dat' / \
        f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}_residuals.csv'
    output_path = config.moments_path / 'moments_dat' / f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}_full.csv'
    output_path_cleaned = config.moments_path / 'moments_dat' / f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}_cleaned.csv'

    # TODO: change clade path
    clade_path = config.intermediate_data_path / 'moments_species_clades_{}.csv'.format(''.join(pops))

    moments_results = moments_utils.load_best_fit_moments_results(inferred_path)

    residual_stats = pd.read_csv(residual_path)
    residual_stats.set_index('species', inplace=True)
    clade_stats = pd.read_csv(clade_path)
    clade_stats.rename(columns={'species_names': 'species'}, inplace=True)
    clade_stats.set_index('species', inplace=True)

    moments_results = rescale_time(moments_results)
    moments_results = moments_results.join(clade_stats).join(residual_stats)

    # now apply a few filters
    print("Print {} passed genome filter".format(moments_results.shape[0]))

    # filter to species with clade structure
    passed_species = moments_results[~moments_results['if_clade']]
    print("Print {} passed clade filter".format(passed_species.shape[0]))

    # filter to species with low uncertainty
    uncert_ratio = passed_species['Tsplit_uncert'] / passed_species['Tsplit']

    # if the uncertainty on the split time is too large, drop the species
    passed = (uncert_ratio < 0.5) & (np.isnan(uncert_ratio)==False)
    passed_species = passed_species[passed]
    print("Print {} passed uncertainty filter".format(sum(passed)))

    resid_threshold = 2
    passed_species = passed_species[passed_species['mean_abs_resid'] < resid_threshold]
    passed_species = passed_species[passed_species['mean_abs_resid_Hadza'] < resid_threshold]
    passed_species = passed_species[passed_species['mean_abs_resid_Tsimane'] < resid_threshold]
    print("Print {} passed residual filter".format(passed_species.shape[0]))

    # rename columns for clarity
    moments_results = moments_results.rename(
        columns={
            'uncert_nu1': 'nu1_std',
            'uncert_nu2': 'nu2_std',
            'uncert_T': 'T_std',
            'uncert_m': 'm_std',
            'uncert_theta': 'theta_std',
            'num_Hadza': 'num_Hadza_MAGs',
            'num_Tsimane': 'num_Tsimane_MAGs',
            'proj_Hadza': 'Hadza_projection',
            'proj_Tsimane': 'Tsimane_projection',
            'syn_genome_length': 'num_total_syn_sites',
            'Tsplit': 'Tsplit (yr)',
            'Tsplit_uncert': 'Tsplit_std (yr)'
        }
    )
    moments_results.to_csv(output_path)
    # now clean up the dataframe for saving
    columns_to_keep = ['nu1', 'nu1_std', 'nu2', 'nu2_std', 'T', 'T_std', 'm', 'm_std', 'theta', 'theta_std', 'num_Hadza_MAGs', 'num_Tsimane_MAGs', 'Hadza_projection', 'Tsimane_projection', 'num_total_syn_sites', 'num_sites_passing_proj', 'Tsplit (yr)', 'Tsplit_std (yr)']
    moments_results = moments_results[columns_to_keep]

    # save the filtered results
    moments_results.loc[passed_species.index].to_csv(output_path_cleaned)
