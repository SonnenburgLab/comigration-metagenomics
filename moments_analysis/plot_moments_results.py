import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import config
from utils import metadata_utils

def main():
    parser = argparse.ArgumentParser(description='Plot moments results.')
    parser.add_argument('--pops', type=str, nargs=2, required=True, help='Pairs of populations to use')
    parser.add_argument('--model_name', type=str, choices=['split_mig', 'split_no_mig'], required=True, help='Model name')
    parser.add_argument('--sfs_batch', type=str, required=True, help='SFS batch to use')
    args = parser.parse_args()

    pops = args.pops
    model_name = args.model_name
    sfs_batch = args.sfs_batch

    popstr = ''.join(pops)

    results_path = config.moments_path / 'moments_dat' / f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}_full.csv'
    moments_results = pd.read_csv(results_path, index_col=0)
    moments_results = moments_results.sort_values(['if_clade', 'mean_abs_resid'], ascending=True)

    moments_figure_path = config.moments_path / 'moments_figures'
    moments_figure_path.mkdir(exist_ok=True)

    databatch = sfs_batch.split('_')[0]
    metadata = metadata_utils.MetadataHelper(databatch)
    print("Total number of species: {}".format(metadata.species.shape[0]))

    print("Print {} passed genome filter".format(moments_results.shape[0]))

    # filter to species with clade structure
    passed_species = moments_results[~moments_results['if_clade']]
    print("Print {} passed clade filter".format(passed_species.shape[0]))

    # filter to species with low uncertainty
    uncert_ratio = passed_species['T_std'] / passed_species['T']

    # if the uncertainty on the split time is too large, drop the species
    passed = (uncert_ratio < 0.5) & (np.isnan(uncert_ratio)==False)
    passed_species = passed_species[passed]
    print("Print {} passed uncertainty filter".format(sum(passed)))

    passed_species_after_uncert = passed_species.copy()
    resid_threshold = 2
    passed_species = passed_species[passed_species['mean_abs_resid'] < resid_threshold]
    passed_species = passed_species[passed_species['mean_abs_resid_Hadza'] < resid_threshold]
    passed_species = passed_species[passed_species['mean_abs_resid_Tsimane'] < resid_threshold]
    print("Print {} passed residual filter".format(passed_species.shape[0]))


    ################################################################
    # first plot residual distribution
    ################################################################
    fig, axes = plt.subplots(1, 2, figsize=(5, 3),)  # Adjust width for better visibility
    plt.subplots_adjust(wspace=0.5)

    # Panel 1: Mean residual size
    sns.boxplot(
        data=moments_results, 
        y='mean_abs_resid', 
        x='if_clade', 
        hue='if_clade', 
        ax=axes[0]
    )
    num_false = moments_results[moments_results['if_clade'] == False].shape[0]
    num_true = moments_results[moments_results['if_clade'] == True].shape[0]
    axes[0].set_xticks([0, 1])
    axes[0].set_xticklabels([f'False\n(n={num_false})', f'True\n(n={num_true})'])
    axes[0].legend().remove()
    axes[0].set_xlabel('')
    # axes[0].set_title('Mean Residual Size', fontsize=10)
    axes[0].set_ylabel('Mean residual size')

    # Panel 2: Split time (log scale)
    sns.boxplot(
        data=moments_results, 
        y='Tsplit (yr)', 
        x='if_clade', 
        hue='if_clade', 
        ax=axes[1]
    )
    axes[1].set_yscale('log')  # Log scale for the y-axis
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels([f'False\n(n={num_false})', f'True\n(n={num_true})'])
    axes[1].legend().remove()
    axes[1].set_xlabel('')
    # axes[1].set_title('Split Time (Years)', fontsize=10)
    axes[1].set_ylabel('Split time (years)')

    axes[0].set_xlabel("If clade structure")
    axes[1].set_xlabel("If clade structure")
    plt.savefig(moments_figure_path / f'{sfs_batch}_residual_size_by_clades_{popstr}_{model_name}.pdf', bbox_inches='tight')
    plt.close()

    ################################################################
    # second plot residual among species with no clades
    ################################################################
    fig, axes = plt.subplots(3, 1, figsize=(6, 6))

    to_plot = passed_species_after_uncert
    to_plot['passed'] = to_plot.index.isin(passed_species.index)

    ax = axes[0]
    sns.barplot(data=to_plot, x='species', y='mean_abs_resid', hue='passed',
                linewidth=1.5, edgecolor=".5", palette=['white', 'tab:blue'], ax=axes[0])
    # ax.axhline(resid_threshold, color='k', linestyle='--')
    ax.set_xticks([])
    ax.set_ylabel('Mean residual size')
    ax.set_xlabel('')

    ax = axes[1]
    sns.barplot(data=to_plot, x='species', y='mean_abs_resid_Tsimane', hue='passed',
                linewidth=1.5, edgecolor=".5", palette=['white', 'tab:blue'], ax=axes[1])
    ax.axhline(resid_threshold, color='k', linestyle='--')
    ax.set_xticks([])
    ax.set_ylabel('Mean residual size\n(Tsimane)')
    ax.set_xlabel('')

    ax = axes[2]
    sns.barplot(data=to_plot, x='species', y='mean_abs_resid_Hadza', hue='passed',
                linewidth=1.5, edgecolor=".5", palette=['white', 'tab:blue'], ax=axes[2])
    ax.axhline(resid_threshold, color='k', linestyle='--')
    # rotate x ticks
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    for ax in axes:
        ax.legend().remove()
    ax.set_ylabel('Mean residual size\n(Hadza)')
    plt.xlabel('')
    plt.savefig(moments_figure_path / f'{sfs_batch}_residual_size_by_species_{popstr}_{model_name}.pdf', bbox_inches='tight')
    plt.close()


    ################################################################
    # third plot split time values
    ################################################################
    plt.figure(figsize=(6, 2))

    # load human parameters
    human_model = pd.read_csv(config.intermediate_data_path / 'MedinaMunoz_model.csv', index_col=0)
    human_model.iloc[:, 1:] = human_model.iloc[:, 1:].astype(float)

    plot_results = passed_species

    # first plot bacterial results
    bacterial_locs = np.array(range(plot_results.shape[0])) + 3
    plt.errorbar(bacterial_locs, plot_results['Tsplit (yr)'] / 1000, yerr=2*plot_results['Tsplit_std (yr)'] / 1000, fmt='o', color='tab:blue')
    bacterial_names = plot_results.index.to_list()

    # next plot human results
    loc = 0
    values = ['TB', 'TN']
    colors = ['tab:blue', 'tab:blue']
    human_locs = []
    for pair in zip(values, colors):
        value = pair[0]
        color = pair[1]
        # plt.errorbar(loc-0.2, human_model.loc[value, 'intronic_value'], yerr=2*human_model.loc[value, 'intronic_SE'], fmt='o', color=color)
        # plt.errorbar(loc+0.2, human_model.loc[value, 'intergenic_value'], yerr=2*human_model.loc[value, 'intergenic_SE'], fmt='o', color=color)
        human_locs.append(loc)
        # plotting only inferred values using synonymous sites
        plt.errorbar(loc, human_model.loc[value, 'syn_value'], yerr=2*human_model.loc[value, 'syn_SE'], fmt='o', color=color)
        plt.axhline(human_model.loc[value, 'syn_value'], color='grey', linestyle='--')
        loc += 1

    plt.xticks(human_locs + list(bacterial_locs), ['Out-of-Africa split', 'Asia-Americas split'] + bacterial_names, 
               rotation=-45, ha='left')
    plt.axvspan(-1, (max(human_locs) + min(bacterial_locs))/2, color='gray', alpha=0.2)

    plt.ylabel('Time (years ago)')
    plt.yscale('log')
    plt.ylim(10, 1e3)
    plt.yticks([10, 100, 1000])
    plt.gca().set_yticklabels(['10,000', '100,000', '1,000,000'])
    plt.xlim(-1, max(bacterial_locs)+1)
    plt.savefig(moments_figure_path / f'{sfs_batch}_split_time_values_{popstr}_{model_name}.pdf', bbox_inches='tight')
    plt.close()

    # now clean up the dataframe for saving
    output_path_cleaned = config.moments_path / 'moments_dat' / f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}_cleaned.csv'
    columns_to_keep = ['nu1', 'nu1_std', 'nu2', 'nu2_std', 'T', 'T_std', 'm', 'm_std', 'theta', 'theta_std', 'num_Hadza_MAGs', 'num_Tsimane_MAGs', 'Hadza_projection', 'Tsimane_projection', 'num_total_syn_sites', 'num_sites_passing_proj', 'Tsplit (yr)', 'Tsplit_std (yr)']
    moments_results = moments_results[columns_to_keep]
    # save the filtered results
    moments_results.loc[passed_species.index].to_csv(output_path_cleaned)

if __name__ == "__main__":
    main()