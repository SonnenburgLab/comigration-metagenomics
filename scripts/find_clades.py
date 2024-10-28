import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from utils import metadata_utils, pairwise_utils
import config

# first genome cutoff
genome_cutoff = 30

metadata = metadata_utils.MetadataHelper(data_batch=config.databatch)
mag_counts = metadata.species[['Hadza', 'Tsimane']]
passed = mag_counts[(mag_counts >= genome_cutoff).all(axis=1)].copy()
passed_species = passed.index.values

pw_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)

clade_species = []
with PdfPages(config.fig_path / 'species_with_clades_div_heatmap.pdf') as pdf:
    for species in passed_species:
        species_helper = pw_helper.get_species_helper(species_name=species)
        clade1, clade2 = species_helper.get_clades(ani_type='ANI', allowed_pops=['Hadza', 'Tsimane'])
        if len(clade1) < 2 or len(clade2) < 2:
            # the minor cluster is not really a clade, just a single MAG that
            # is distant from the rest
            continue

        # for printing ANI distributions
        # clade1_ani = species_helper.get_ANI_dist_by_mags(clade1)
        # clade2_ani = species_helper.get_ANI_dist_by_mags(clade2)
        # clade12_ani = species_helper.get_ANI_dist_between_mags(clade1, clade2)
        # print(np.mean(clade1_ani), np.mean(clade2_ani), np.mean(clade12_ani))

        g = species_helper.visualize_clades(ani_type='ANI', allowed_pops=['Hadza', 'Tsimane'])
        pdf.savefig()
        plt.close()
        clade_species.append(species)

passed['if_clade'] = passed.index.isin(clade_species)
passed.to_csv(config.intermediate_data_path / 'moments_species.csv')