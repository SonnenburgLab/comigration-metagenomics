import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import scipy
from scipy.cluster.hierarchy import linkage, fcluster
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch

import config
from utils import pairwise_utils, metadata_utils

# full_res_annotated = pd.read_csv(config.run_path / '240625_all_runs_all_species_annotated.csv')
full_res_annotated = pd.read_csv(config.run_path / f'{config.databatch}_annotated.csv')
metadata = metadata_utils.MetadataHelper(config.databatch)
all_assignments = []
fig_dir = Path('figs/240730_divergence_clustermap.pdf')
with PdfPages(fig_dir) as pdf:
    for species_name, species_res in full_res_annotated.groupby('species'):
        print(species_name)
        # first, prepare pairwise table
        # Use combine_first to fill NaN values in the pivot_table with those in transpose_table
        symmetric_table = pairwise_utils.long_form_to_symmetric(species_res, row_name='genome1', col_name='genome2', val_name='div')

        # second, remove closely related pairs
        mean_div = species_res['div'].mean()
        dedup_genomes = pairwise_utils.cluster_close_pairs(symmetric_table, cut_dist=mean_div / 10)
        symmetric_table = symmetric_table.loc[dedup_genomes, dedup_genomes]

        # third, make dendrogram with the remaining data
        condensed_y = scipy.spatial.distance.squareform(symmetric_table)
        Z = linkage(condensed_y, method='average')  # You can change 'average' to other methods like 'single', 'complete', etc.

        # optional: keeping only the major clade by clustering in to 2 clusters
        # max_clusters = 2
        # clusters = fcluster(Z, max_clusters, criterion='maxclust')
        # cluster_assignments = pd.Series(clusters, index=symmetric_table.index)
        # all_assignments.append(cluster_assignments)
        # # unused function to color the heatmap by cluster
        # palette = sns.color_palette("pastel", max_clusters) 
        # color_dict = {i+1: color for i, color in enumerate(palette)}
        # cluster_colors = [color_dict[cluster] for cluster in clusters]

        # finally, plot the heatmap
        # coloring by population
        pop_assignments = [metadata.get_mag_pop(x) for x in symmetric_table.index]
        all_pops = metadata.get_all_pops()
        palette = sns.color_palette("pastel", len(all_pops))  # You can choose different palettes
        color_dict = {pop: color for pop, color in zip(all_pops, palette)}
        # make paleofeces stand out more
        color_dict['Paleofeces'] = 'black'
        row_colors = [color_dict[pop] for pop in pop_assignments]

        g = sns.clustermap(symmetric_table, row_linkage=Z, col_linkage=Z, row_colors=row_colors, cmap='viridis', vmin=0, vmax=0.08, cbar_kws=dict(orientation='horizontal'))

        # position the legend for population coloring
        x0, y0, w, h = g.cbar_pos
        g.ax_cbar.set_position([x0, 0.95, g.ax_row_dendrogram.get_position().width * 0.9, 0.02])

        handles = [Patch(facecolor=color_dict[name]) for name in all_pops]
        plt.legend(handles, all_pops, title='Population', 
                   bbox_to_anchor=(1.05, -0.05), bbox_transform=g.ax_heatmap.transAxes, loc='upper left')

        plt.title(species_name)
        pdf.savefig(bbox_inches='tight')
        plt.close()

# all_assignments = pd.concat(all_assignments)
# all_assignments.to_csv(config.run_path / '240625_divergence_clustermap_assignments.csv')