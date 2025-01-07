import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
import time

from identical_tract_analysis.simulate_mutation_accumulation import get_passed_species_full_within
from utils import pairwise_utils
import config

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
        print(f"Warning: L99_target ({L99_target}) is less than the minimum L99 in the simulation ({sim_L99s.min()})")
        return sim_years[-1], None
    if L99_target > sim_L99s[0]:
        print(f"Warning: L99_target ({L99_target}) is greater than the maximum L99 in the simulation ({sim_L99s.max()})")
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

    num_bootstrap = 100
    pair_threshold = 200
    all_focal_pops = [['Hadza', 'Tsimane'], ['China', 'HMP'], ['HMP', 'MetaHIT']]
    sim_path = config.intermediate_data_path / 'mutation_accumulation'
    sim_years = np.logspace(3, 5.5, 20)

    res_data = []
    for species in species_list:
        species_helper = pairwise_helper.get_species_helper(species)
        species_sim_data = pd.read_csv(sim_path / f'{species}_full_mutation_accumulation.tsv', sep='\t')

        for focal_pops in all_focal_pops:
            between_l99, between_size = species_helper.get_between_pop_L99(focal_pops[0], focal_pops[1])
            print(species, focal_pops, between_l99, between_size)
            if between_size < pair_threshold:
                continue
            inferred_years = np.zeros(num_bootstrap)
            for i in range(num_bootstrap):
                bs_l99s = get_simulation_bootstrap_L99(species_sim_data, between_size, len(sim_years))
                # smooth the L99s to be isotonic
                bs_l99s_smoothed = iso_reg.fit_transform(sim_years, bs_l99s)
                inferred_year, func = find_L99_year(sim_years, bs_l99s_smoothed, between_l99)
                inferred_years[i] = inferred_year
            mean_year = np.mean(inferred_years)
            ci = np.percentile(inferred_years, [2.5, 97.5])
            res_data.append([species, focal_pops[0], focal_pops[1], between_l99, between_size, mean_year, ci[0], ci[1]])

    res_df = pd.DataFrame(res_data, columns=['species', 'pop1', 'pop2', 'L99', 'num_pairs', 'mean_year', 'ci_low', 'ci_high'])