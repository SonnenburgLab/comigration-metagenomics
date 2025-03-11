"""
Take the inferred model parameters and compute the residuals of model fits
"""
import pandas as pd
import numpy as np
import argparse
import moments
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import config
from utils import moments_utils

def compute_moments_residuals(pops=['Hadza', 'Tsimane'], model_name='split_mig', sfs_batch=config.sfs_batch):
    # loading raw inference results
    inferred_path = config.moments_path / 'moments_dat' / \
        f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}.csv'

    moments_results = moments_utils.load_best_fit_moments_results(inferred_path)

    # compute the residual of model predictions
    residuals = {}
    if model_name == 'split_mig':
        model_func = moments.Demographics2D.split_mig
    elif model_name == 'split_no_mig':
        model_func = moments_utils.split_no_mig

    residual_stats = pd.DataFrame(index=moments_results.index, columns=['mean_abs_resid', 'mean_abs_resid_Tsimane', 'mean_abs_resid_Hadza'])
    for i, species in enumerate(moments_results.index):
        print(f"Computing residuals for {species} ({i+1}/{len(moments_results)})")
        # load data, prepare a model SFS, and compute residuals
        data, _ = moments_utils.load_SFS_projection(species, focal_pops=pops, sfs_folder=config.sfs_path / sfs_batch)
        try:
            model = moments_utils.prep_model(moments_results.loc[species], 
                                            model_func, data, model_name=model_name)
        except IndexError:
            print("Failed to prep model for {}".format(species))
            continue

        resid = moments.Inference.linear_Poisson_residual(model, data, mask=2)
        residuals[species] = resid

        # save mean resid size
        flatresid = np.array(np.compress(np.logical_not(resid.mask.ravel()), resid.ravel()))
        residual_stats.loc[species, 'mean_abs_resid'] = np.abs(flatresid).mean()

        # also compute the marginal residuals
        f_h = data.marginalize([1])
        f_t = data.marginalize([0])
        model_h = model.marginalize([1])
        model_t = model.marginalize([0])
        resid_t = moments.Inference.linear_Poisson_residual(model_t, f_t)
        resid_h = moments.Inference.linear_Poisson_residual(model_h, f_h)
        # compute mean for SNVs with > 1 alt allele
        mean_resid_t = np.abs(resid_t[2:]).mean()
        mean_resid_h = np.abs(resid_h[2:]).mean()
        residual_stats.loc[species, 'mean_abs_resid_Tsimane'] = mean_resid_t
        residual_stats.loc[species, 'mean_abs_resid_Hadza'] = mean_resid_h

    # save flat residual values
    # currently not used in the analysis

    # resid_dfs = []
    # for species, resid in residuals.items():
    #     flatresid = np.array(np.compress(np.logical_not(resid.mask.ravel()), resid.ravel()))
    #     resid_df = pd.DataFrame(flatresid)
    #     resid_df.columns = ['residual']
    #     resid_df['species'] = species
    #     resid_dfs.append(resid_df)
    # full_resid_df = pd.concat(resid_dfs)
    return residual_stats

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute residuals of model fits.')
    parser.add_argument('--pops', type=str, nargs=2, required=True, help='Pairs of populations to use')
    parser.add_argument('--model_name', type=str, choices=['split_mig', 'split_no_mig'], required=True, help='Model name')
    parser.add_argument('--sfs_batch', type=str, required=True, help='SFS batch to use')
    args = parser.parse_args()

    pops = args.pops
    model_name = args.model_name
    sfs_batch = args.sfs_batch

    residual_stats = compute_moments_residuals(pops=pops, model_name=model_name, sfs_batch=sfs_batch)
    residual_stats.to_csv(config.moments_path / 'moments_dat' / f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}_residuals.csv')
    print("Done!")