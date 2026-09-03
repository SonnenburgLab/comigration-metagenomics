import time
import os

from utils import pairwise_utils
import config

# cphmm is now an installed package (`pip install -e` the close_pair_hmm repo); see README.
import cphmm.prior
import cphmm.infer_pipelines as infer_pipelines
import tsimane_datahelper

pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)

species_list = pairwise_helper.get_species_list()
result_path = config.cphmm_res_path / 'results'
result_path.mkdir(parents=True, exist_ok=True)

for species in species_list:
    print("Processing species {} at {}".format(species, time.ctime()))
    if not os.path.exists(cphmm.prior.get_prior_filename(species, prior_path=config.cphmm_prior_path)):
        print("Skipping species {} because of lack of prior".format(species))
        continue
    summary_file = os.path.join(result_path, species + '__summary.csv')
    transfers_file = os.path.join(result_path, species + '__transfers.csv')
    if os.path.exists(summary_file):
        print("Skipping species {} because summary file exists".format(species))
        continue

    species_dat = pairwise_helper.drep_summary[pairwise_helper.drep_summary['species']==species]
    dh = tsimane_datahelper.DataHelper_Hadza_Tsimane(species=species, drep_summary=species_dat)
    # use our own generated priors, not the package's bundled set
    dh.hmm_prior_path = config.cphmm_prior_path

    infer_summary, transfer_summary = infer_pipelines.infer_pairs(dh, dh.get_close_pairs(perc_id_threshold=0.5))

    infer_summary.to_csv(summary_file)
    transfer_summary.to_csv(transfers_file)
    print("Done processing species {} at {}".format(species, time.ctime()))

print("Done at {}".format(time.ctime()))
