import time
import os
import pickle

from utils import pairwise_utils
import config

# cphmm is now an installed package (`pip install -e` the close_pair_hmm repo); see README.
import cphmm.prior
import tsimane_datahelper

pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)

species_list = pairwise_helper.get_species_list()
print("Total species: ", len(species_list))

error_species = {}
for species in species_list:
    print("Processing species {} at {}".format(species, time.ctime()))
    if os.path.exists(cphmm.prior.get_prior_filename(species, prior_path=config.cphmm_prior_path)):
        print("Skipping species {} as prior already exists".format(species))
        continue
    species_dat = pairwise_helper.drep_summary[pairwise_helper.drep_summary['species']==species]
    dh = tsimane_datahelper.DataHelper_Hadza_Tsimane(species=species, drep_summary=species_dat)

    try:
        local_divs, genome_divs = cphmm.prior.sample_blocks(dh)
    except ValueError as e:
        error_species[species] = e
        print("Skipping species {} due to ValueError".format(species))
        continue
    divs, counts = cphmm.prior.compute_div_histogram(local_divs, genome_divs, separate_clades=False)
    cphmm.prior.save_prior(divs, counts, dh.species, prior_path=config.cphmm_prior_path)

pickle.dump(error_species, open('error_species.pkl', 'wb'))
