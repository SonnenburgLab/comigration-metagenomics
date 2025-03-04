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
    mut_per_year = mut_rate * day_per_year * gen_per_day
    return div / mut_per_year

infer_summary = pd.read_csv(config.cphmm_res_path / '241022_inference_summary_full.tsv', sep='\t')
infer_summary.set_index(['genome1', 'genome2'], inplace=True)

# combine the perc id data
pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)
drep_summary = pairwise_helper.hgt_summary[pairwise_helper.hgt_summary['species'].isin(infer_summary['species'].unique())]
drep_summary = drep_summary.set_index(['genome1', 'genome2'])

infer_summary['ani'] = drep_summary.loc[infer_summary.index, 'ani']
infer_summary['perc_id'] = drep_summary.loc[infer_summary.index, 'perc_id']
infer_summary['div_years_binned'] = div_to_years(infer_summary['est_div'])
infer_summary['div_years'] = div_to_years(infer_summary['naive_div'])


# fit trend per species
species_trend = pd.DataFrame(columns=['species', 'b', 'year_at_10perc', 'num_pairs', 'genome_len'])
with PdfPages('figs/241022_percid_trend/241022_perc_id_species_trend.pdf') as pdf:
    for species, dat in infer_summary.groupby('species'):
        # filter 1: number of pairs
        if len(dat) < 20:
            continue

        # filter 2: not too many zero divergence pairs
        if (dat['div_years_binned'] == 0).sum() > 0.3 * len(dat):
            continue

        year, b = fit_perc_id_trend(dat['perc_id'], dat['div_years_binned'])
        species_trend.loc[species] = [species, b, year, len(dat), dat['genome_len'].mean()]

        plt.subplots(figsize=(4, 3))
        x = dat['div_years_binned']
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

species_trend.to_csv('figs/241022_percid_trend/241022_perc_id_species_trend.tsv', sep='\t')

# hist of years
plt.figure()
sns.histplot(species_trend['year_at_10perc'], bins=20)
plt.axvline(species_trend['year_at_10perc'].median(), color='grey', linestyle='--', 
            label="Median: {}".format(species_trend['year_at_10perc'].median()))
plt.xlabel('Divergence years at 10% identical genes')
plt.ylabel('Number of species')
plt.legend()
plt.savefig('figs/241022_percid_trend/241022_perc_id_years_hist.pdf')


# next, plot the aggregated trend
species_to_include = species_trend['species'].unique()
infer_summary_filtered = infer_summary[infer_summary['species'].isin(species_to_include)]

year, b = fit_perc_id_trend(infer_summary_filtered['perc_id'], infer_summary_filtered['div_years_binned'])

plt.subplots(figsize=(6, 4))
x = infer_summary_filtered['div_years_binned']
y = infer_summary_filtered['perc_id']

yplot = np.linspace(0, 1)
plt.plot(prec_id_trend_func(yplot, b), yplot, '-', color='red', 
        label="Average trend\nacross {} pairs".format(len(x)))
plt.scatter(year, 0.1, label='10% ~ {} years'.format(int(year)),
            color='red', zorder=10)
plt.plot([0, year], [0.1, 0.1], color='grey', linestyle='--')
plt.plot([year, year], [0, 0.1], color='grey', linestyle='--')

binx = np.linspace(0, 10000, 100)
biny = np.linspace(0, 1, 100)

sns.histplot(x=x, y=y, bins=(binx, biny), cbar=True, norm=LogNorm(vmin=10, vmax=1e3), vmin=None, vmax=None,
             cmap='viridis')

# plt.axhspan(0, 0.5, color='grey', alpha=0.1)
plt.fill_between([0, 10000], [0.5, 0.5], color='grey', alpha=0.2, hatch='//', label='Data not used')

plt.xlim(0, 10000)
plt.ylim(0, 1)
plt.xlabel('Clonal divergence (years)')
plt.ylabel('Percent identical genes')
plt.legend(loc='lower right')

plt.savefig('figs/241022_percid_trend/241022_perc_id_trend.pdf', bbox_inches='tight')