import numpy as np
import pandas as pd
import scipy

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

"""
For loading data
"""

def proj_func(counts):
    # for determining the number of samples to project to
    return int(np.median(counts) * 0.95)

def load_SFS_projection(species, proj_func=proj_func, sfs_folder=config.sfs_path / config.databatch, focal_pops=['Hadza', 'Tsimane']):
    """
    Loading the SFS and projecting to a smaller size using a projection function
    """
    sfs_file = f'{sfs_folder}/{species}.snps.txt'
    dd = dadi.Misc.make_data_dict(sfs_file)
    hz_counts = [sum(dd[x]['calls'][focal_pops[0]]) for x in dd]
    ts_counts = [sum(dd[x]['calls'][focal_pops[1]]) for x in dd]

    proj = [proj_func(hz_counts), proj_func(ts_counts)]
    data = moments.Spectrum.from_data_dict(dd, pop_ids=focal_pops, projections=proj, polarized=False)
    return data

def load_SFS(species, metadata, proj_ratio=0.9, sfs_folder=config.sfs_path / config.databatch, focal_pops=['Hadza', 'Tsimane']):
    """
    Loading the SFS while saving some stats about how many sites are included
    """
    mag_count_df = metadata.get_mag_counts(species)
    # choosing a projection size; corresponding to a prevalence cutoff
    proj = mag_count_df.loc[focal_pops].values.astype(int) * proj_ratio
    proj = proj.astype(int)
    return _load_SFS(species, proj, sfs_folder=sfs_folder, focal_pops=focal_pops)

def _load_SFS(species, proj, sfs_folder=config.sfs_path / config.databatch, focal_pops=['Hadza', 'Tsimane']):
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


"""
Functions for analyzing fitted models
"""

# sfs_folder = '/Users/Device6/Documents/Research/bgoodlab/microbiome_codiv/comigration_metagenomics/dat/240801_dadi_dict'
def prep_model(result_row, model_func, data, model_name='split_mig'):
    """
    Take the result row from the csv of inferred model parameters
    and simulate SFS under the demographic model

    Needs data to infer the optimal theta for comparison
    """
    if model_name == 'split_mig':
        params = result_row[['nu1', 'nu2', 'T', 'm']]
    elif model_name == 'split_no_mig':
        params = result_row[['nu1', 'nu2', 'T']]
    theta = result_row['theta']
    model = model_func(params, data.sample_sizes)
    data.mask[:, 1] = True
    data.mask[1, :] = True
    theta = moments.Inference.optimal_sfs_scaling(model, data)
    model = model * theta
    model = model.fold()
    data.mask[:, 1] = False
    data.mask[1, :] = False
    return model

def prep_model_unscaled(result_row, model_func, sample_sizes=(100, 100)):
    params = result_row[['nu1', 'nu2', 'T', 'm']]
    theta = result_row['theta']
    model = model_func(params, sample_sizes)
    model = model.fold()
    return model

def prep_model_neutral(fs):
    """
    Prepare the neutral single pop SFS that matches the sample size and theta
    of observed SFS
    """
    model = moments.Demographics1D.snm(fs.sample_sizes)
    fs.mask[1] = True
    opt_theta = moments.Inference.optimal_sfs_scaling(model, fs)
    model = model * opt_theta
    model = model.fold()
    fs.mask[1] = False
    return model

"""
Various ways of parsing the saved SFS data
"""

def load_raw_SFS_proj(species, proj, sfs_folder=config.sfs_path / config.databatch, focal_pops=['Hadza', 'Tsimane']):
    sfs_file = f'{sfs_folder}/{species}.snps.txt'
    dd = dadi.Misc.make_data_dict(sfs_file)

    data = moments.Spectrum.from_data_dict(dd, pop_ids=focal_pops, projections=proj, polarized=False)
    return data

def load_sfs_counts(species, sfs_folder=config.sfs_path / config.databatch, pops=['Hadza', 'Tsimane']):
    # loading the raw counts from dadi dicts
    sfs_file = f'{sfs_folder}/{species}.snps.txt'
    dd = dadi.Misc.make_data_dict(sfs_file)

    all_counts = []
    snv_ids = []
    for snv in dd:
        snv_ids.append(snv)
        snv_dat = []
        for pop in pops:
            pop_dat = dd[snv]['calls'][pop]
            snv_dat.append(pop_dat[0])
            snv_dat.append(pop_dat[1])
        all_counts.append(snv_dat)

    all_counts = np.array(all_counts)

    sfs_dat = pd.DataFrame(all_counts, columns=['Hadza_ref', 'Hadza_alt', 'Tsimane_ref', 'Tsimane_alt'], index=snv_ids)

    sfs_dat['Hadza_tot'] = sfs_dat['Hadza_alt'] + sfs_dat['Hadza_ref']
    sfs_dat['Tsimane_tot'] = sfs_dat['Tsimane_alt'] + sfs_dat['Tsimane_ref']
    sfs_dat['Total_ref'] = sfs_dat['Hadza_ref'] + sfs_dat['Tsimane_ref']
    sfs_dat['Total_alt'] = sfs_dat['Hadza_alt'] + sfs_dat['Tsimane_alt']

    # polarize to major allele
    mask = sfs_dat['Total_ref']<sfs_dat['Total_alt']
    rows = sfs_dat[mask].copy()
    sfs_dat.loc[mask, 'Hadza_ref'] = rows['Hadza_alt']
    sfs_dat.loc[mask, 'Hadza_alt'] = rows['Hadza_ref']
    sfs_dat.loc[mask, 'Tsimane_ref'] = rows['Tsimane_alt']
    sfs_dat.loc[mask, 'Tsimane_alt'] = rows['Tsimane_ref']
    sfs_dat.loc[mask, 'Total_ref'] = rows['Total_alt']
    sfs_dat.loc[mask, 'Total_alt'] = rows['Total_ref']

    sfs_dat['Is_exclusive'] = (sfs_dat['Hadza_alt']==0) | (sfs_dat['Tsimane_alt']==0)

    # drop monomorphic sites
    sfs_dat = sfs_dat[sfs_dat['Total_alt']>0]
    return sfs_dat


"""
Functions for computing likelihoods of private SNVs
"""

def well_mixed_likelihood(k1, k_tot, n1, n2):
    """
    The likelihood for how allele counts distribute in two populations
    Assuming that the two populations are well mixed, so described by
    a binomial distribution.

    :param k1: number of observed alleles in population 1
    :param k_tot: total number of observed alleles
    :param n_h: total number of haplotypes in population 1
    :param n_t: total number of haplotypes in population 2
    """
    n_tot = n1 + n2
    p = n1 / n_tot
    return scipy.stats.binom.pmf(k1, k_tot, p)


def well_mixed_exclusive_snv_likelihood(k_tot, n1, n2):
    """
    The likelihood that all observed alleles are from a single population
    Assuming that the two populations are well mixed, so described by
    a binomial distribution.

    :param k_tot: total number of observed alleles
    :param n1: total number of haplotypes in population 1
    :param n2: total number of haplotypes in population 2
    """
    return well_mixed_likelihood(k_tot, k_tot, n1, n2) + well_mixed_likelihood(0, k_tot, n1, n2)


def compute_well_mixed_likelihoods(sfs_dat):
    sfs_dat['Likelihood'] = well_mixed_likelihood(sfs_dat['Hadza_alt'], sfs_dat['Total_alt'], 
                                sfs_dat['Hadza_ref']+sfs_dat['Hadza_alt'], sfs_dat['Tsimane_ref']+sfs_dat['Tsimane_alt'])

    sfs_dat['Exclusive_snv_likelihood'] = well_mixed_exclusive_snv_likelihood(sfs_dat['Total_alt'],
                                            sfs_dat['Hadza_alt']+sfs_dat['Hadza_ref'], sfs_dat['Tsimane_alt'] + sfs_dat['Tsimane_ref'])
    
    return sfs_dat


def compute_exclusive_likelihood_model(sfs_dat, model):
    model_probs = pd.DataFrame(index=sfs_dat.index, columns=['Model_likelihood', 'Model_exclusive_likelihood'], dtype=float)
    for sizes, grouped in sfs_dat.groupby(['Hadza_tot', 'Tsimane_tot']):
        proj_model = moments.Spectrum.project(model, sizes)
        for snv_id, row in grouped.iterrows():
            # sizes = row['Hadza_tot'], row['Tsimane_tot']
            k1 = row['Hadza_alt']
            k2 = row['Tsimane_alt']
            k_tot = row['Total_alt']
            exclu_prob = 0
            if k_tot < sizes[0]:
                exclu_prob += proj_model[k_tot, 0]
            if k_tot < sizes[1]:
                exclu_prob += proj_model[0, k_tot]

            # sum over all possible configurations that k1+k2 = k_tot for computing the total weight of observing k_tot
            total_prob = 0
            for i in range(0, k_tot+1):
                if i < (sizes[0]+1) and (k_tot-i < sizes[1]+1) and k_tot-i >= 0:
                    total_prob += proj_model[i, k_tot-i]

            model_probs.loc[snv_id, 'Model_exclusive_likelihood'] = exclu_prob / total_prob
            model_probs.loc[snv_id, 'Model_likelihood'] = proj_model[k1, k2] / total_prob
    return model_probs