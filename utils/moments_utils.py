import numpy as np

import demes
import moments
import dadi

import config

"""
A few custom demographic models for moments.
"""
def split_no_mig(params, ns):
    nu1, nu2, T_div = params
    return moments.Demographics2D.split_mig((nu1, nu2, T_div, 0), ns)

def split_mig_iso(params, ns):
    """
    Split with migration model with full isolation between populations after a while.
    """
    nu1, nu2, m, T1, T2 = params

    # f for the equilibrium ancestral population
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)

    # The divergence
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])

    fs.integrate([nu1, nu2], T1, m=np.array([[0, m], [m, 0]]))
    fs.integrate([nu1, nu2], T2, m=np.array([[0, 0], [0, 0]]))

    return fs

def growth_mixed(params, ns):
    """
    Improved exp growth model that takes in rate directly
    """
    rate, T = params

    nu_func = lambda t: [np.exp(rate * t)]
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate(nu_func, T, 0.01)

    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    return fs

def growth_split_mig(growth_param, params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must have length 2")
    
    rate, T = growth_param 
    Tsplit, m = params

    nuB = 1
    nu_anc = lambda t: [nuB * np.exp(rate * t / T)]
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate(nu_anc, T - Tsplit, dt_fac=0.01)
    # we split the population
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    # size at split
    nu0 = nu_anc(T - Tsplit)[0]
    nu_func = lambda t: 2 * [nu0 * np.exp(rate * t / T)]
    # continue to grow after split
    fs.integrate(nu_func, Tsplit, m=np.array([[0, m], [m, 0]]))
    fs.pop_ids = pop_ids
    return fs

def growth_split(growth_param, params, ns):
    return growth_split_mig(growth_param, (params[0], 0), ns)

"""
Wrapper functions for inferences
"""

def fit_model(data, model_func, p_guess, lower_bound, upper_bound, mask_singletons=True):

    if mask_singletons:
        data.mask[1, :] = True
        data.mask[:, 1] = True
    
    p_guess = moments.Misc.perturb_params(
        p_guess, lower_bound=lower_bound, upper_bound=upper_bound)

    opt_params = moments.Inference.optimize_log_fmin(
        p_guess, data, model_func,
        lower_bound=lower_bound, upper_bound=upper_bound,
        verbose=20) # report every 20 iterations

    refit_theta = moments.Inference.optimal_sfs_scaling(
        model_func(opt_params, data.sample_sizes), data)

    uncerts = moments.Godambe.FIM_uncert(
        model_func, opt_params, data)

    # print_params = params + [input_theta]
    print_opt = np.concatenate((opt_params, [refit_theta]))
    ll = moments.Inference.ll(model_func(opt_params, data.sample_sizes)*refit_theta, data)
    return opt_params, refit_theta, uncerts, ll


def load_SFS(species, metadata, proj_ratio=0.9, focal_pops=['Hadza', 'Tsimane']):
    data_batch = metadata.data_batch
    mag_count_df = metadata.get_mag_counts(species)
    # choosing a projection size; corresponding to a prevalence cutoff
    proj = mag_count_df.loc[focal_pops].values.astype(int) * proj_ratio
    proj = proj.astype(int)
    return _load_SFS(species, proj, data_batch, focal_pops)

def _load_SFS(species, proj, data_batch, focal_pops=['Hadza', 'Tsimane']):
    sfs_folder = config.sfs_path / data_batch

    sfs_file = sfs_folder / f'{species}.snps.txt'

    dd = dadi.Misc.make_data_dict(sfs_file)
    try:
        data = moments.Spectrum.from_data_dict(dd, pop_ids=focal_pops, projections = proj, polarized = False)
    except AttributeError:
        # because of stupid python version problem, need this way to catch the data warning in moments
        # this error is because no SNV is found after projection
        return None

    num_total_snvs = len(dd)
    num_focal_snvs = count_segragating_snvs_in_pops(dd, focal_pops)
    num_proj_snvs = data.S()
    stats = {
        'num_total_snvs': num_total_snvs,
        'num_focal_snvs': num_focal_snvs,
        'num_proj_snvs': num_proj_snvs
    }
    for i, pop in enumerate(focal_pops):
        stats[f'proj_{pop}'] = proj[i]
    return data, stats


def count_segragating_snvs_in_pops(data_dict, pops):
    all_is_seg = []
    for entry in data_dict:
        is_seg = []
        for pop in pops:
            is_seg.append(data_dict[entry]['calls'][pop][1]>0)
        all_is_seg.append(is_seg)
    all_is_seg = np.array(all_is_seg)
    num_seg_in_focal = (all_is_seg.any(axis=1)).sum()
    return num_seg_in_focal

# def get_SFS_stats(species, data_batch):
#     sfs_file = config.sfs_path / data_batch / f'{species}.snps.txt'
#     dd = dadi.Misc.make_data_dict(sfs_file)

def rescale_time(T, theta, num_sites, mu=4.08e-10, gen_per_day=1):
    """
    In the model, theta=4 * N_ref * mu, while the unit of T is 2 * N_ref
    """
    N_anc = theta / num_sites / mu / 2
    t_gen = T * N_anc
    t_years = t_gen / (365 * gen_per_day)
    return t_years