import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

import config
from utils import pairwise_utils


def prec_id_trend_func(y, b):
    # y is the percentage of identical genes
    # assuming a linear decrease over time
    # this returns the predicted clonal divergence age of the pair
    return b * (1 - y)

def fit_perc_id_trend(perc_id_vals, clonal_years):
    popt, pcov = curve_fit(prec_id_trend_func, perc_id_vals, clonal_years)
    year_at_10perc = popt[0] * 0.9
    return year_at_10perc, popt[0]

def div_to_years(div, gen_per_day=config.gen_per_day, mut_rate=config.mut_rate, day_per_year=config.day_per_year):
    # factor of 2 because of pairwise comparison
    mut_per_year = 2 * mut_rate * day_per_year * gen_per_day
    return div / mut_per_year

def load_infer_summary(pairwise_helper=None):
    # loading cphmm inference results
    infer_summary = pd.read_csv(config.cphmm_res_path / '241022_inference_summary_full.tsv', sep='\t')
    infer_summary.set_index(['genome1', 'genome2'], inplace=True)
    print(f"Old total number of pairs: {len(infer_summary)}")

    # combine the perc id data
    if pairwise_helper is None:
        pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)
    drep_summary = pairwise_helper.drep_summary[pairwise_helper.drep_summary['species'].isin(infer_summary['species'].unique())]
    drep_summary = drep_summary.set_index(['genome1', 'genome2'])

    # find intersection of index
    shared_index = infer_summary.index.intersection(drep_summary.index)
    print(f"Remaining total number of pairs: {len(shared_index)}")
    infer_summary = infer_summary.loc[shared_index]

    infer_summary['ani'] = drep_summary.loc[infer_summary.index, 'ani']
    infer_summary['perc_id'] = drep_summary.loc[infer_summary.index, 'perc_id']
    infer_summary['div_years'] = div_to_years(infer_summary['naive_div'])
    # a more conservative estimate of clonal divergence, which leads to higher recombination rate estimates
    infer_summary['div_years_binned'] = div_to_years(infer_summary['est_div'])
    return infer_summary


if __name__ == '__main__':
    infer_summary = load_infer_summary()
    infer_summary.drop(columns=['est_div'], inplace=True)
    infer_summary.to_csv(config.supp_table_path / 'cphmm_infer_summary.tsv', sep='\t')

    ################################################################
    # fit trend per species
    ################################################################
    fig_path = config.project_path / 'cphmm_analysis' / 'close_pair_figs'
    species_trend = pd.DataFrame(columns=['species', 'b', 'year_at_10perc', 'num_pairs', 'genome_len'])
    fit_val = 'div_years'
    with PdfPages(fig_path / 'perc_id_species_trend.pdf') as pdf:
        for species, dat in infer_summary.groupby('species'):
            # filter 1: number of pairs
            if len(dat) < 20:
                continue

            # filter 2: not too many zero divergence pairs
            if (dat[fit_val] == 0).sum() > 0.3 * len(dat):
                continue

            year, b = fit_perc_id_trend(dat['perc_id'], dat[fit_val])
            species_trend.loc[species] = [species, b, year, len(dat), dat['genome_len'].mean()]

            plt.subplots(figsize=(4, 3))
            x = dat[fit_val]
            y = dat['perc_id']

            plt.scatter(x, y, alpha=0.5)
            yplot = np.linspace(0, 1)
            plt.plot(prec_id_trend_func(yplot, b), yplot, '-', color='grey', 
                    label='10% ~ {} years'.format(int(year)))

            plt.scatter(year, 0.1, color='grey')
            plt.axvline(x=year, color='grey', linestyle='--')
            plt.axhline(y=0.1, color='grey', linestyle='--')

            plt.ylim([0, 1])
            plt.xlim([0, 15000])

            plt.xlabel('Clonal divergence (years)')
            plt.ylabel('Percent identical genes')
            plt.title(species)
            plt.legend()
            pdf.savefig(bbox_inches='tight')
            plt.close()

    species_trend.to_csv(config.project_path / 'cphmm_analysis' / 'perc_id_species_trend.tsv', sep='\t')


    ################################################################
    # plot examples
    ################################################################

    species_list = ['Acetatifactor_acetigignens', 'Akkermansia_muciniphila', 
                    'Bulleidia_sp905194705', 'Faecalibacterium_prausnitzii_M', 
                    'Prevotella_sp900548745', 'Ruminococcus_E_bromii_B']

    fig, axes = plt.subplots(3, 2, figsize=(6, 6))

    plt.subplots_adjust(hspace=0.5, wspace=0.3)

    for i, species in enumerate(species_list):
        ax = axes[i//2, i%2]

        dat = infer_summary[infer_summary['species'] == species]
        year, b = fit_perc_id_trend(dat['perc_id'], dat[fit_val])
        species_trend.loc[species] = [species, b, year, len(dat), dat['genome_len'].mean()]

        y = dat[fit_val]
        x = dat['perc_id']

        ax.scatter(1-x, y, alpha=0.5, s=3)
        xplot = np.linspace(0, 1)
        ax.plot(1-xplot, prec_id_trend_func(xplot, b), '-', color='grey',)

        ax.scatter(1 - 0.1, year, s=20, color='grey', label='10% ~ {} years'.format(int(year)))
        ax.axhline(y=year, color='grey', linestyle='--')
        ax.axvline(x=1 - 0.1, color='grey', linestyle='--')
        ax.fill_between([0.5, 1],[10000, 10000],  color='grey', alpha=0.1, hatch='//')
        ax.set_xticklabels(['100', '80', '60', '40', '20', '0'])
        ax.legend(loc='upper left')
        ax.set_title(species)

        ax.set_xlim([0, 1])
        ax.set_ylim([0, 5000])

    for i in range(2):
        axes[i, 0].set_xticklabels([])
        axes[i, 1].set_xticklabels([])
    for ax in axes[:, 1]:
        ax.set_yticklabels([])

    # shared x label
    fig.text(0.5, 0.02, 'Fraction of identical genes (%)', ha='center', fontsize=10)
    fig.text(0.01, 0.5, 'Clonal divergence (years)', va='center', rotation='vertical', fontsize=10)

    fig_path = config.supp_fig_path / 'perc_id_trend_examples.pdf'
    fig.savefig(fig_path, bbox_inches='tight')

    ################################################################
    # plot inferred years histogram
    ################################################################
    plt.figure(figsize=(4, 3))
    sns.histplot(species_trend['year_at_10perc'], bins=30)
    median_val = species_trend['year_at_10perc'].median()
    plt.axvline(median_val, color='tab:orange', linestyle='--', label='median={0:.0f} years'.format(median_val))
    plt.xlabel('Clonal divergence at \n10% identiacal fraction (years)')
    plt.legend()
    plt.savefig(config.supp_fig_path / 'perc_id_trend_hist.pdf', bbox_inches='tight')


    ################################################################
    # next, plot the aggregated trend
    ################################################################

    species_to_include = species_trend['species'].unique()
    infer_summary_filtered = infer_summary[infer_summary['species'].isin(species_to_include)]

    year, b = fit_perc_id_trend(infer_summary_filtered['perc_id'], infer_summary_filtered[fit_val])
    print("Aggregated constant of proportionality (1/years)", 1/b)

    # set small font size
    plt.rcParams.update({'font.size': 8})
    plt.subplots(figsize=(5, 3))
    y = infer_summary_filtered[fit_val]
    x = infer_summary_filtered['perc_id']

    xplot = np.linspace(0, 1)
    plt.plot(1 - xplot, prec_id_trend_func(xplot, b), '-', color='tab:orange', 
            label="Average trend\nacross {} pairs".format(len(x)))
    plt.scatter(1 - 0.1, year, label='10% ~ {} years'.format(int(year)),
                color='tab:orange', zorder=10)
    plt.plot([1 - 0.1, 1 - 0.1], [year, 0], color='grey', linestyle='--')
    plt.plot([0, 1 - 0.1], [year, year], color='grey', linestyle='--')

    biny = np.linspace(0, 5000, 100)
    binx = np.linspace(0, 1, 100)

    sns.histplot(x=1-x, y=y, bins=(binx, biny), cbar=True, norm=LogNorm(vmin=10, vmax=1e3), vmin=None, vmax=None,
                cmap='Blues')

    # plt.axhspan(0, 0.5, color='grey', alpha=0.1)
    plt.fill_between([0.5, 1],[10000, 10000],  color='grey', alpha=0.1, hatch='//', label='Detection limit')
    plt.gca().set_xticklabels(['100', '80', '60', '40', '20', '0'])

    plt.ylim(0, 5000)
    plt.xlim(0, 1)
    plt.ylabel('Clonal divergence (years)')
    plt.xlabel('Fraction of identical genes (%)')
    plt.legend(loc='upper right')

    plt.savefig(config.supp_fig_path / 'perc_id_trend_full.pdf', bbox_inches='tight')