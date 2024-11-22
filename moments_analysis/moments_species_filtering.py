import pandas as pd
import numpy as np

import moments

import config
from utils import moments_utils

def filter_species_for_pops(pops=['Hadza', 'Tsimane'], model_name='split_mig'):
    popstr = ''.join(pops)
    # first load the clade annotation
    moments_species = pd.read_csv(config.intermediate_data_path / 'moments_species_clades_{}.csv'.format(popstr), index_col=0)

    # loading demographic inference results
    moment_results = pd.read_csv(config.project_path / 'moments_out' / 'unclustered' / '240714__{}__{}__{}.csv'.format(model_name, pops[0], pops[1]))
    moment_results.set_index('species', inplace=True)

    # use mutation rate to translate scaled time to years
    moment_results['Tsplit'] = moments_utils.rescale_time(moment_results['T'], moment_results['theta'], moment_results['num_sites_passing_proj'])
    moment_results['Tsplit_uncert'] = moments_utils.rescale_time(moment_results['uncert_T'], moment_results['theta'], moment_results['num_sites_passing_proj'])

    # compute the residual of model predictions
    residuals = {}
    if model_name == 'split_mig':
        model_func = moments.Demographics2D.split_mig
    elif model_name == 'split_no_mig':
        model_func = moments_utils.split_no_mig

    for i, species in enumerate(moments_species.index):
        proj = moment_results.loc[species, ['proj_'+pop for pop in pops]].astype(int).to_list()
        
        # data, _ = moments_utils._load_SFS(species, proj, config.databatch)
        data = moments_utils.load_SFS_projection(species, focal_pops=pops)
        try:
            model = moments_utils.prep_model(moment_results.loc[species], 
                                            model_func, data, model_name=model_name)
        except IndexError:
            print("Failed to prep model for {}".format(species))
            continue

        resid = moments.Inference.linear_Poisson_residual(model, data, mask=1)
        residuals[species] = resid

    # save mean residual size (unsigned)
    residual_stats = pd.DataFrame(index=moments_species.index, columns=['mean_abs_resid'])
    resid_dfs = []
    for species, resid in residuals.items():
        flatresid = np.array(np.compress(np.logical_not(resid.mask.ravel()), resid.ravel()))
        resid_df = pd.DataFrame(flatresid)
        resid_df.columns = ['residual']
        resid_df['species'] = species
        resid_dfs.append(resid_df)

        residual_stats.loc[species, 'mean_abs_resid'] = np.abs(flatresid).mean()

    # might be useful if we want to plot the distribution of residuals
    full_resid_df = pd.concat(resid_dfs)

    # combine with other species info
    moments_species = moments_species.join(residual_stats['mean_abs_resid']).sort_values(
        ['if_clade', 'mean_abs_resid'], ascending=True)
    # moments_species = moments_species.join(moment_results[['Tsplit', 'Tsplit_uncert']])
    moments_species = moments_species.join(moment_results)

    moments_species.to_csv(config.project_path / 'moments_analysis' / 'moments_results_{}_{}.csv'.format(popstr, model_name))


if __name__ == '__main__':
    # filter_species_for_pops(pops=['Hadza', 'Tsimane'], model_name='split_no_mig')
    filter_species_for_pops(pops=['MetaHIT', 'HMP'], model_name='split_no_mig')
    print("Done!")