import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from sklearn.isotonic import IsotonicRegression
import time

from identical_tract_analysis.simulate_mutation_accumulation import get_passed_species_full_within
from utils import pairwise_utils
import config

def join_dfs_by_genome_pair(df1, df2):
    # Helper function to join two dataframes on the 'genome1' and 'genome2' columns
    # Create new columns where genome1 and genome2 are sorted alphabetically
    # join type is inner, meaning rows not in df1 are dropped
    df1['genome_pair'] = df1.apply(lambda row: tuple(sorted([row['genome1'], row['genome2']])), axis=1)
    df2['genome_pair'] = df2.apply(lambda row: tuple(sorted([row['genome1'], row['genome2']])), axis=1)

    # Drop the original genome1 and genome2 (optional)
    df1 = df1.drop(columns=['genome1', 'genome2'])
    df2 = df2.drop(columns=['genome1', 'genome2'])

    # Now join the two dataframes on the 'genome_pair' column
    merged_df = pd.merge(df1, df2, on='genome_pair', how='inner')

    # Split the 'genome_pair' tuple back into two columns for MultiIndex
    merged_df[['genome1', 'genome2']] = pd.DataFrame(merged_df['genome_pair'].tolist(), index=merged_df.index)

    # Drop the 'genome_pair' column (optional)
    merged_df = merged_df.drop(columns=['genome_pair'])

    # Set the new 'genome1' and 'genome2' columns as a MultiIndex
    merged_df.set_index(['genome1', 'genome2'], inplace=True)
    return merged_df

def get_simulation_bootstrap_L99(sim_data, num_to_sample, num_years):
    L99_bootstrap = np.zeros(num_years)
    for j in range(num_years):
        resampled = sim_data['sim_max_run_{}'.format(j)].sample(n=num_to_sample, replace=True)
        L99_bootstrap[j] = resampled.quantile(0.99)
    return L99_bootstrap

def find_x_for_y(target_y, interp_func, x_bounds, x0, x1):
    """
    Find the x value at which the interpolated function equals target_y
    """
    def func_to_solve(x_val):
        return interp_func(x_val) - target_y

    # Solve for x in the range of x_bounds
    result = root_scalar(func_to_solve, bracket=x_bounds, x0=x0, x1=x1)

    if result.converged:
        return result.root
    else:
        raise ValueError("Root finding did not converge.")

def find_L99_year(sim_years, sim_L99s, L99_target, log_scale=True):
    """
    Take a L99 vs year curve and find the year at which the L99 is closest to the target L99

    If log_scale is True, the linear interpolation is done in log space
    """
    if L99_target < sim_L99s[-1]:
        # print(f"Warning: L99_target ({L99_target}) is less than the minimum L99 in the simulation ({sim_L99s.min()})")
        return sim_years[-1], None
    if L99_target > sim_L99s[0]:
        # print(f"Warning: L99_target ({L99_target}) is greater than the maximum L99 in the simulation ({sim_L99s.max()})")
        return sim_years[0], None
    # Step 1: Interpolate a function (linear or spline interpolation)
    if log_scale:
        interp_func = interp1d(np.log10(sim_years), np.log10(sim_L99s), kind='linear', fill_value='extrapolate')
        interp_func_ret = lambda x: 10**interp_func(np.log10(x))
        L99_target = np.log10(L99_target)
        x_bounds = (3, 5.5)
        x0 = 3
        x1 = 5.5
    else:
        interp_func = interp1d(sim_years, sim_L99s, kind='linear')
        interp_func_ret = interp_func
        x_bounds = (1e3, 3e5)
        x0 = 1e3
        x1 = 3e5

    # Provide bounds for the search (min and max of your x range)
    try:
        x_result = find_x_for_y(L99_target, interp_func, x_bounds, x0=x0, x1=x1)
        if log_scale:
            x_result = 10**x_result
        # print(f"x value for y = {L99_target} is approximately {x_result}")
        return x_result, interp_func_ret
    except ValueError as e:
        print(e)
        return None, interp_func_ret

if __name__ == '__main__':
    species_list = get_passed_species_full_within()
    pairwise_helper = pairwise_utils.PairwiseHelper(databatch=config.databatch)

    num_bootstrap = 1000
    # if a population comparison has fewer than this number of pairs, skip it
    pair_threshold = 200

    all_focal_pops = [['Hadza', 'Tsimane'], ['China', 'HMP'], ['HMP', 'MetaHIT']]
    sim_path = config.intermediate_data_path / 'mutation_accumulation'

    # simulated range of years in simulate_mutation_accumulation.py
    sim_years = np.logspace(3, 5.5, 20)

    # Fit isotonic regression
    iso_reg = IsotonicRegression(increasing=False)

    res_data = []
    for species in species_list:
        species_helper = pairwise_helper.get_species_helper(species)
        species_sim_data = pd.read_csv(sim_path / f'{species}_full_mutation_accumulation.tsv', sep='\t')
        species_sim_data = join_dfs_by_genome_pair(
            species_sim_data, species_helper.run_summary[['genome1', 'genome2', 'max_run']].copy())

        for focal_pops in all_focal_pops:
            between_l99, between_size = species_helper.get_between_pop_L99(focal_pops[0], focal_pops[1])
            if between_size < pair_threshold:
                continue
            print(species, focal_pops, between_l99, between_size)
            # first compute within population L99s and CI
            within_l99_bs = pairwise_utils.compute_L99_resampled(species_sim_data,
                                val_name='max_run', n_bootstrap=num_bootstrap,
                                n_resample=between_size)
            within_l99_bs = np.array(within_l99_bs)
            mean_within_l99 = within_l99_bs.mean()
            std_within_l99 = within_l99_bs.std()

            # next compute p val vs the shortest simulation time
            sim_0_l99_bs = pairwise_utils.compute_L99_resampled(species_sim_data,
                                val_name='sim_max_run_0', n_bootstrap=num_bootstrap,
                                n_resample=between_size)
            sim_0_l99_bs = np.array(sim_0_l99_bs)
            # test if shows significant isolation longer than the shortest simulation time
            p_val = np.sum(sim_0_l99_bs < between_l99) / num_bootstrap
            if p_val > 0.05:
                print("No significant isolation > 1kyr")

            # next, infer isolation years from bootstrap L99 year curves
            inferred_years = np.zeros(num_bootstrap)
            for i in range(num_bootstrap):
                bs_l99s = get_simulation_bootstrap_L99(species_sim_data, between_size, len(sim_years))
                # smooth the L99s to be isotonic
                bs_l99s_smoothed = iso_reg.fit_transform(sim_years, bs_l99s)
                inferred_year, func = find_L99_year(sim_years, bs_l99s_smoothed, between_l99)
                inferred_years[i] = inferred_year
            mean_year = np.mean(inferred_years)
            ci = np.percentile(inferred_years, [2.5, 97.5])
            res_data.append([species, focal_pops[0], focal_pops[1], between_l99, between_size, mean_year, ci[0], ci[1], p_val, mean_within_l99, std_within_l99])

    res_df = pd.DataFrame(res_data, columns=['species', 'pop1', 'pop2', 'between_L99', 'num_pairs', 'mean_year', 'ci_low', 'ci_high', 'p_val', 'mean_within_L99', 'std_within_L99'])
    res_df.to_csv(config.intermediate_data_path / 'mutation_accumulation' / 'between_pop_inferred_years.tsv', sep='\t', index=False)