# CP-HMM analysis

This folder contains the code to run CP-HMM inference of recombination events on
closely related MAGs. The main inference package was developed for [Liu & Good (2024)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002472) and can be found in this [repository](https://github.com/zhiru-liu/close_pair_hmm).

- `tsimane_datahelper.py` contains a wrapper class of the `snv_utils.PairwiseSNVHelper` class. It implements a few functions necessary for interfacing with the CP-HMM package, such as returning the synonymous SNV profile between a pair of MAGs.
- CP-HMM inference requires a prior distribution of transfer divergences, which is computed by sampling genome blocks from random pairs of MAGs. `prior_prep.py` runs this step for the Tsimane dataset.
- `infer_all.py` runs the inference on all pair of closely related MAGs in the Tsimane dataset and caches the results to `config.cphmm_res_path / 'results'`. For each bacterial species, the inference produces a pairwise summary file containing the total number of recombination events and the inferred clonal diveregence for each pair of MAGs, as well as a transfer file containing the location of each recombination event.
    - This step can be time consuming if there's a large number of closely related MAGs. Please consider adapting the script to parallelize the inference per species on a cluster.
    - `concat_results.py` simply concatenates the results of the inference into a single file.
- `fit_divergence_years.py` takes in the pairwise summary file and fits a linear relationship between divergence years and the percentage of identical genes. The results are used to estimate the divergence time between a pair of almost fully recombined MAGs (i.e. <10% identical genes.)
    - This script also generates a few summary plots, such as a histogram of the inferred divergence times.