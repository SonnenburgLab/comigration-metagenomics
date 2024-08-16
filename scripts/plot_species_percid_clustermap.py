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

def plot_clustermap(species_list, fig_dir, pairwise_helper=None):
    if pairwise_helper is None:
        pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)
    metadata = metadata_utils.MetadataHelper(config.databatch)
    with PdfPages(fig_dir) as pdf:
        for species in species_list:
            print(species)
            species_helper = pairwise_helper.get_species_helper(species)
            symmetric_table = 1 - species_helper.get_percid_mat()
            Z = pairwise_utils.linkage_clustering(symmetric_table)
            clusters = species_helper.get_clonal_clusters()

            cluster_sizes = clusters.groupby('cluster').size()
            big_clusters = cluster_sizes[cluster_sizes > 1]
            colors = sns.color_palette('tab20', len(big_clusters))
            cluster_colors = dict(zip(big_clusters.index, colors))
            clusters['cluster_color'] = clusters['cluster'].map(lambda x: cluster_colors.get(x, 'white'))

            # finally, plot the heatmap
            # coloring by population
            pop_assignments = [metadata.get_mag_pop(x) for x in symmetric_table.index]
            all_pops = metadata.get_all_pops()
            palette = sns.color_palette("pastel", len(all_pops))  # You can choose different palettes
            color_dict = {pop: color for pop, color in zip(all_pops, palette)}
            # make paleofeces stand out more
            color_dict['Paleofeces'] = 'black'
            row_colors = [color_dict[pop] for pop in pop_assignments]

            g = sns.clustermap(1 - symmetric_table, row_linkage=Z, col_linkage=Z, row_colors=row_colors, col_colors=clusters['cluster_color'], \
                cmap='viridis', vmin=0, vmax=1, cbar_kws=dict(orientation='horizontal'), figsize=(10, 10))

            # position the legend for population coloring
            x0, y0, w, h = g.cbar_pos
            g.ax_cbar.set_position([x0, 0.95, g.ax_row_dendrogram.get_position().width * 0.9, 0.02])

            handles = [Patch(facecolor=color_dict[name]) for name in all_pops]
            plt.legend(handles, all_pops, title='Population', 
                    bbox_to_anchor=(1.05, -0.05), bbox_transform=g.ax_heatmap.transAxes, loc='upper left')

            plt.title(species)
            pdf.savefig(bbox_inches='tight')
            plt.close()

if __name__ == '__main__':
    # full_res_annotated = pd.read_csv(config.run_path / f'{config.databatch}_annotated.csv')
    # species_list = full_res_annotated['species'].unique()
    # fig_dir = Path('figs/240730_divergence_clustermap.pdf')

    fig_dir = Path('figs/240812_percid_clustermap.pdf')
    pairwise_helper = pairwise_utils.PairwiseHelper(config.databatch)
    species_list = pairwise_helper.get_species_list()
    plot_clustermap(species_list, fig_dir, pairwise_helper=pairwise_helper)