# Moments analysis steps

This folder contains the scripts for Moments demographic analysis.

## Prepare data

- Compute SFS using the SNV catalogs. Saved in dadi dict format.
Use script `save_sfs_dict.py`
    - If this is the first time loading the SNV catalogs, `snv_utils.SNVHelper` will take some time to compute the 4D sites.

## Snakemake Pipeline

The moments analysis has been automated using Snakemake. The pipeline includes the following steps:

1. Fit Moments models using `fit_moments.py`.
2. Compute residual statistics using `compute_residuals.py`.
3. Postprocess results using `postprocess_results.py`.
4. Plot the full SFS using `plot_moments_SFS_full.py`.
5. Generate diagnostic plots using `plot_moments_results.py`.

The pipeline is defined in the `Snakefile` and can be executed by running:
```sh
snakemake -j <number_of_jobs>
```