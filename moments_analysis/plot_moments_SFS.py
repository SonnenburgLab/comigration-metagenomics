import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import moments

from utils import moments_utils
from moments_analysis.plot_moments_SFS_full import plot_one_species
import config

sfs_batch = config.sfs_batch
model_name = 'split_mig'
pops = ['Hadza', 'Tsimane']
moment_res_path = config.moments_path / 'moments_dat' / f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}_cleaned.csv'
# moment_res_path = config.supp_table_path / 'supp_moments_results_HadzaTsimane_split_mig.tsv'
moments_results = pd.read_csv(moment_res_path, index_col=0)

font_size = 3
mpl.rcParams['axes.titlesize'] = font_size
mpl.rcParams['axes.labelsize'] = font_size
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size
mpl.rcParams['legend.fontsize'] = font_size


passed_species = moments_results.index.to_list()

####################################################################################################
# Plot SFS and residuals for the first 8 species
####################################################################################################
species_to_plot = passed_species[:8]
fig, axes = plt.subplots(len(species_to_plot), 6, figsize=(6.5, 8))
plt.subplots_adjust(hspace=0.7, wspace=0.5)

for i, species in enumerate(species_to_plot):
    plot_one_species(axes[i, :], moments_results, species)
axes[0, 0].legend()

fig.savefig(config.moments_path / 'moments_figures' / 'moments_sfs_1.pdf', dpi=300, bbox_inches='tight')
plt.close()

####################################################################################################
# Plot SFS and residuals for the remaining species
####################################################################################################
species_to_plot = passed_species[8:]
fig, axes = plt.subplots(len(species_to_plot), 6, figsize=(6.5, len(species_to_plot)))
plt.subplots_adjust(hspace=0.7, wspace=0.5)

for i, species in enumerate(species_to_plot):
    plot_one_species(axes[i, :], moments_results, species)
axes[0, 0].legend()

# fig.savefig(config.supp_fig_path / 'moments_sfs_2.pdf', dpi=300, bbox_inches='tight')
fig.savefig(config.moments_path / 'moments_figures' / 'moments_sfs_2.pdf', dpi=300, bbox_inches='tight')
plt.close()