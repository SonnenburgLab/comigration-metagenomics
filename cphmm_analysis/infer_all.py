import sys
import time
import os

from utils import pairwise_utils
import config

# manually adding the cphmm path
sys.path.append('/Users/Device6/Documents/Research/bgoodlab/close_pair_hmm')
import cphmm.prior
import tsimane_datahelper
import infer_pipelines

pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)

species_list = pairwise_helper.get_species_list()
result_path = config.cphmm_res_path / 'results'

for species in species_list:
    print("Processing species {} at {}".format(species, time.ctime()))
    if not os.path.exists(cphmm.prior.get_prior_filename(species)):
        print("Skipping species {} because of lack of prior".format(species))
        continue
    summary_file = os.path.join(result_path, species + '__summary.csv')
    transfers_file = os.path.join(result_path, species + '__transfers.csv')
    if os.path.exists(summary_file):
        print("Skipping species {} because summary file exists".format(species))
        continue

    species_dat = pairwise_helper.drep_summary[pairwise_helper.drep_summary['species']==species]
    dh = tsimane_datahelper.DataHelper_Hadza_Tsimane(species=species, drep_summary=species_dat)

    infer_summary, transfer_summary = infer_pipelines.infer_pairs(dh, dh.get_close_pairs(perc_id_threshold=0.5))

    infer_summary.to_csv(summary_file)
    transfer_summary.to_csv(transfers_file)
    print("Done processing species {} at {}".format(species, time.ctime()))

print("Done at {}".format(time.ctime()))