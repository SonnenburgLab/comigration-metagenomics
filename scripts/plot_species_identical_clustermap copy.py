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
hgt_res = pd.read_csv('/Users/Device6/Documents/Research/bgoodlab/microbiome_codiv/dRep_analysis/Global_HGT_v7.1.csv')
hgt_res['species'] = hgt_res['name'].str.replace('_v7', '')
hgt_res['genome1'] = hgt_res['querry'].str.replace('\.fa', '')
hgt_res['genome2'] = hgt_res['reference'].str.replace('\.fa', '')

metadata = metadata_utils.MetadataHelper(config.databatch)
all_assignments = []
fig_dir = Path('figs/240730_perc_id_clustermap.pdf')
with PdfPages(fig_dir) as pdf:
    for species_name, species_res in hgt_res.groupby('species'):
        print(species_name)
        # first, prepare pairwise table
        # Use combine_first to fill NaN values in the pivot_table with those in transpose_table
        try:
            symmetric_table = pairwise_utils.long_form_to_symmetric(species_res, row_name='genome1', col_name='genome2', val_name='perc_id')
        except ValueError:
            print(f'No HGT results for {species_name}')
            continue
        np.fill_diagonal(symmetric_table.values, 1)
        symmetric_table = 1 - symmetric_table

        condensed_y = scipy.spatial.distance.squareform(symmetric_table)
        full_Z = linkage(condensed_y, method='average')  # You can change 'average' to other methods like 'single', 'complete', etc.

        # now cut tree
        passed_mags, clusters = pairwise_utils.cluster_close_pairs(symmetric_table, 0.8, return_clusters=True)
        passed_pd_mat = symmetric_table.loc[passed_mags, passed_mags]

        condensed_y = scipy.spatial.distance.squareform(passed_pd_mat)
        Z = linkage(condensed_y, method='average')  # You can change 'average' to other methods like 'single', 'complete', etc.

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

        g = sns.clustermap(1 - symmetric_table, row_linkage=full_Z, col_linkage=full_Z, row_colors=row_colors, col_colors=clusters['cluster_color'], \
            cmap='viridis', vmin=0, vmax=1, cbar_kws=dict(orientation='horizontal'), figsize=(10, 10))

        # position the legend for population coloring
        x0, y0, w, h = g.cbar_pos
        g.ax_cbar.set_position([x0, 0.95, g.ax_row_dendrogram.get_position().width * 0.9, 0.02])

        handles = [Patch(facecolor=color_dict[name]) for name in all_pops]
        plt.legend(handles, all_pops, title='Population', 
                   bbox_to_anchor=(1.05, -0.05), bbox_transform=g.ax_heatmap.transAxes, loc='upper left')

        plt.title(species_name)
        pdf.savefig(bbox_inches='tight')
        plt.close()

        clusters['cluster'].to_csv(f'dat/clusters/240731_{species_name}_clusters.csv')