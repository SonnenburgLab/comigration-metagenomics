import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import config

def div_to_years(div, gen_per_day=1, mut_rate=4.08e-10):
    mut_per_year = 2 * mut_rate * 365 * gen_per_day
    return div / mut_per_year

# combine the perc id data
from utils import pairwise_utils
pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)
drep_summary = pairwise_helper.hgt_summary[pairwise_helper.hgt_summary['species'].isin(infer_summary['species'].unique())]
drep_summary = drep_summary.set_index(['genome1', 'genome2'])
mean_syn_div = pairwise_helper.pairwise_summary.groupby('species')['div'].mean()
mean_ani = pairwise_helper.hgt_summary.groupby('species')['ani'].mean()
pairwise_div_pd = pd.DataFrame({'mean_syn_div': mean_syn_div, 'mean_ani': mean_ani})
pairwise_div_pd['ani_diff'] = 1 - pairwise_div_pd['mean_ani']

species_trend = pd.read_csv(config.intermediate_data_path / 'cphmm_species_years.csv')
species_summary = species_trend.join(pairwise_div_pd, on='species')
species_summary['year_at_Tc'] = div_to_years(species_summary['mean_syn_div'])
species_summary['year_at_Tm'] = species_summary['year_at_10perc'] / 0.9

# load Liu Good 2024 data
plosbio_Tm = pd.read_csv('/Volumes/Botein/LiuGood2024_files/figs_cleanup/supp_table/TcTm_estimation.csv')

fig, ax = plt.subplots(figsize=(4, 3))
bins = np.logspace(0, 2, 20)

sns.histplot(species_summary['year_at_Tc'] / species_summary['year_at_Tm'], bins=bins, stat='density', ax=ax)
sns.histplot(plosbio_Tm['Tc_years'] / plosbio_Tm['Tm_years'], bins=bins, stat='density', ax=ax)
ax.set_xscale('log')
n_species = len(species_summary)
n_plosbio = len(plosbio_Tm)
ax.set_xlabel('$T_{mrca} / T_{mosaic}$')
ax.set_xlim([1, 300])
ax.legend(['This study (n={})'.format(n_species), 'Liu Good 2024 (n={})'.format(n_plosbio)])