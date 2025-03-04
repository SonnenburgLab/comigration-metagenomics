# Moments analysis steps

This folder contains the scripts for Moments demographic analysis.

## Steps

- Compute SFS using the SNV catalogs. Saved in dadi dict format.
Use script `save_sfs_dict.py`
    - If this is the first time loading the SNV catalogs, `snv_utils.SNVHelper` will take some time to compute the 4D sites.
- Fit Moments models using `fit_moments.py`
    - If haven't computed the site coverage yet, run `data_processing/compute_and_save_coverage_matrix.py` first

### Filtering and quality control
- Compute residual statistics using `compute_residuals.py`
- Identify potential clade structure using `find_clades.py`
    - Requires pairwise ANI data between MAGs
- Filter species based on a few filters using `postprocess_results.py`
<!-- - Run `moments_species_filtering.py` to compute residual statistics -->

### Plotting SI figures
- Run `plot_moments_results.py` to plot diagnostic plots + count number of species filtered at each step
- After finding the final set of passed species shown in main text, compute private SNV statistics with `plot_private_SNVs.py`