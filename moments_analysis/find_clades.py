"""
Identifying species that may have clade structure in the Hadza and Tsimane 
populations. We use hierarchical clustering to identify the two clusters that
are most distinct from each other. We then check if the minor cluster has at
least 2 MAGs. If so, we consider this species to potentially have a clade 
structure. Note this is a very simple heuristic and may not be perfect, but
is a conservative way to filter out clade species for moments analysis.

This script will also generate a heatmap for each species that fall into
this category. The heatmap will show the pairwise ANI distances between MAGs
within the species and help double check if the clade structure looks "real".
"""
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from utils import metadata_utils, pairwise_utils
import config

def find_clade_for_pops(genome_cutoff=30, pops=['Hadza', 'Tsimane']):
    metadata = metadata_utils.MetadataHelper(data_batch=config.databatch)
    mag_counts = metadata.species[pops]
    passed = mag_counts[(mag_counts >= genome_cutoff).all(axis=1)].copy()
    passed_species = passed.index.values

    pw_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)

    clade_species = []
    popstr = ''.join(pops)
    with PdfPages(config.project_path / 'moments_analysis' / 'species_with_clades_div_heatmap_{}.pdf'.format(popstr)) as pdf:
        for species in passed_species:
            species_helper = pw_helper.get_species_helper(species_name=species)
            clade1, clade2 = species_helper.get_clades(ani_type='ANI', allowed_pops=pops)
            if len(clade1) < 2 or len(clade2) < 2:
                # the minor cluster is not really a clade, just a single MAG that
                # is distant from the rest
                continue

            # for printing ANI distributions
            # clade1_ani = species_helper.get_ANI_dist_by_mags(clade1)
            # clade2_ani = species_helper.get_ANI_dist_by_mags(clade2)
            # clade12_ani = species_helper.get_ANI_dist_between_mags(clade1, clade2)
            # print(np.mean(clade1_ani), np.mean(clade2_ani), np.mean(clade12_ani))

            g = species_helper.visualize_clades(ani_type='ANI', allowed_pops=pops)
            pdf.savefig()
            plt.close()
            clade_species.append(species)

    passed['if_clade'] = passed.index.isin(clade_species)
    passed.to_csv(config.intermediate_data_path / 'moments_species_clades_{}.csv'.format(popstr))


if __name__ == '__main__':
    find_clade_for_pops()
    find_clade_for_pops(pops=['MetaHIT', 'HMP'])
    print("Done!")