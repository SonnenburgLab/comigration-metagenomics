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


full_res_annotated = pd.read_csv(config.run_path / '240627_all_runs_all_species_annotated.csv')
metadata = metadata_utils.MetadataHelper('240626')
all_assignments = []
fig_dir = Path('figs/240627_pop_HGT.pdf')

with PdfPages(fig_dir) as pdf:
    for species_name, species_res in full_res_annotated.groupby('species'):
        print(species_name)
        species_res = full_res_annotated[full_res_annotated['species'] == species_name]
        average_div = species_res['div'].mean()
        species_res = pairwise_utils.dedup_pairwise_df(species_res, cut_dist=average_div / config.close_pair_div_ratio)

        dedup_genomes = np.unique(species_res[['genome1', 'genome2']])

        pop_comp_mat = pairwise_utils.compute_pop_comparisons(species_res, pairwise_utils.HGT_score)
        pair_count_mat = pairwise_utils.compute_pop_comparisons(species_res, len)

        mag_counts = metadata.get_mags_pop_count(dedup_genomes)
        pop_name_map = {x: x + f' ({mag_counts[x]})' for x in mag_counts.to_dict()}
        # filter out populations with less than 5 genomes
        mag_counts = mag_counts[mag_counts > 5]
        pop_comp_mat = pop_comp_mat.loc[mag_counts.index, mag_counts.index]
        pair_count_mat = pair_count_mat.loc[mag_counts.index, mag_counts.index]
        if pop_comp_mat.shape[0] < 2:
            continue

        # plt.figure(figsize=(3, 3))
        g = sns.clustermap(pop_comp_mat.fillna(0), mask=pair_count_mat<100, cmap='coolwarm', vmin=-0.5, vmax=1.5, annot=True,
                        cbar_kws=dict(orientation='horizontal'))
        g.ax_heatmap.set_xticklabels([pop_name_map[label.get_text()] for label in g.ax_heatmap.get_xticklabels()])
        g.ax_heatmap.set_yticklabels([pop_name_map[label.get_text()] for label in g.ax_heatmap.get_yticklabels()])
        x0, y0, w, h = g.cbar_pos
        g.ax_cbar.set_position([x0, 0.85, g.ax_row_dendrogram.get_position().width * 0.9, 0.02])
        g.ax_cbar.set_title(f'{species_name}\n\n\nHGT score')

        pdf.savefig(bbox_inches='tight')
        plt.close()
