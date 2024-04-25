import pandas as pd
import numpy as np
import datetime
import dadi

def write_dadi_input(output_file, snps, pops=['Hadza', 'Tsimane']):
    """
    Write the snps to a file in the format required by dadi
    # https://dadi.readthedocs.io/en/latest/user-guide/importing-data/

    snps is a list of dataframes, each containing the snps for a population in the same order as pops

    Each dataframe should have the following multi-index, which specify a unique snv:
    contig, position, ref, alt

    The values of the dataframe should be 0, 1, or 255, where 0 is the reference allele, 1 is the alternate allele, and 255 is missing data.
    """
    # can have arbitrary number of populations for now
    header_items = ['Ref', 'Out', 'Allele1']
    for pop in pops:
        header_items.append(pop)
    header_items.append('Allele2')
    for pop in pops:
        header_items.append(pop)
    header_items = header_items + ['Contig', 'Position', '\n']
    
    # print time
    print('Writing to file:', output_file)
    with open(output_file, 'w') as snp_file:
        snp_file.write('\t'.join(header_items))

        for ind, hz_row in snps[0].iterrows():
            # iterate over snp rows
            # first collect snp information
            ref_string = '-'+ind[2]+'-'
            out_string = '---'
            a1 = ind[2]
            a2 = ind[3]
            pos = ind[1]
            contig = ind[0]

            refs = []
            alts = []
            # now write over populations
            for pop_snp in snps:
                pop_row = pop_snp.loc[ind, :]
                refs.append(np.sum(pop_row==0))
                alts.append(np.sum(pop_row==1))
            
            line_items = [ref_string, out_string, a1] + refs + [a2] + alts + [contig, pos, '\n']
            line_items = [str(x) for x in line_items]
            snp_file.write('\t'.join(line_items))
    print('Finished at:', datetime.datetime.now())


"""
Parameter fitting functions
"""

def exp_single(params, ns, pts):
    """
    Single population with exponential growth.
    nu: Ratio of contemporary to ancestral population size.
    T: Time in the past at which growth began (in units of 2*Na generations)
    """
    nu, T = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    nu_func = lambda t: np.exp(np.log(nu) * t/T)
    phi = dadi.Integration.one_pop(phi, xx, T, nu_func)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

def no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    """
    nu1, nu2, T = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def no_mig_anc(params, ns, pts):
    """
    Split into two populations, no migration.
    Ancestral population is under exponential growth.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    nuA: Size of ancestral population
    TA: Time in the past at which ancestral population growth began
    """
    nu1, nu2, T, nuA, TA = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # ancestral population growth for TA + T before now
    nu_func = lambda t: np.exp(np.log(nuA) * t/TA)
    phi = dadi.Integration.one_pop(phi, xx, TA, nu_func)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # split population evolve for T gen
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def no_mig_anc_two_epoch(params, ns, pts):
    """
    Split into two populations, no migration.
    Ancestral population is under exponential growth.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    nuA: Size of ancestral population
    TA: Time in the past at which ancestral population growth began
    """
    nu1, nu2, T, nuA, TA = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    
    phi = dadi.Integration.one_pop(phi, xx, TA, nuA)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # split population evolve for T gen
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def no_mig_anc_TS_exp(params, ns, pts):
    """
    Split into two populations, no migration.
    Ancestral population is under exponential growth.
    Then Tsimane population is under exponential growth.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    nuA: Size of ancestral population
    TA: Time in the past at which ancestral population growth began
    nuTS: Size of Tsimane population after exponential growth
    """
    nu1, nu2, T, nuA, TA, nuTS = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # ancestral population growth for TA + T before now
    nu_func = lambda t: np.exp(np.log(nuA) * t/TA)
    phi = dadi.Integration.one_pop(phi, xx, TA, nu_func)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nuTS_func = lambda t: np.exp(np.log(nuTS/nu2) * t/T)
    # split population evolve for T gen
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nuTS_func, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def mig(params, ns, pts):
    """
    Split into two populations, with migration.
    """
    nu1, nu2, T, m = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def no_mig_exp(params, ns, pts):
    """
    Split into two populations, no migration.
    Each population undergoes exponential growth after the split.
    """
    nu1, nu2, nu1_fi, nu2_fi, T = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1 * np.exp(np.log(nu1_fi/nu1) * t/T)
    nu2_func = lambda t: nu2 * np.exp(np.log(nu2_fi/nu2) * t/T)

    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def mig_two_epoch(params, ns, pts):
    """
    Split into two populations, with migration.
    Each population undergoes exponential growth after the split.
    """
    nu1, nu2, nu1_fi, nu2_fi, T, T1, T2, m= params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1 if t < T * T1 else nu1_fi
    nu2_func = lambda t: nu2 if t < T * T2 else nu2_fi

    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=m, m21=m)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def mig_two_epoch_varying_m(params, ns, pts):
    """
    Split into two populations, with migration.
    Each population undergoes two epoch growth after the split.
    Migration stops after Tiso * T generations.
    """
    nu1, nu2, nu1_fi, nu2_fi, T, T1, T2, m, m2, Tiso = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1 if t < T * T1 else nu1_fi
    nu2_func = lambda t: nu2 if t < T * T2 else nu2_fi
    mu_func = lambda t: m if t < T * Tiso else m2

    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=mu_func, m21=mu_func)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def mig_anc_two_epoch(params, ns, pts):
    """
    Split into two populations, no migration.
    Ancestral population is under exponential growth.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    nuA: Size of ancestral population
    TA: Time in the past at which ancestral population growth began
    mu: Migration rate
    """
    nu1, nu2, T, nuA, TA, mu = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    
    phi = dadi.Integration.one_pop(phi, xx, TA, nuA)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # split population evolve for T gen
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=mu, m21=mu)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def mig_exp(params, ns, pts):
    """
    Split into two populations, with migration.
    Each population undergoes exponential growth after the split.
    """
    nu1, nu2, nu1_fi, nu2_fi, T, m = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1 * np.exp(np.log(nu1_fi/nu1) * t/T)
    nu2_func = lambda t: nu2 * np.exp(np.log(nu2_fi/nu2) * t/T)

    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=m, m21=m)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def mig_iso_exp(params, ns, pts):
    """
    Split into two populations, with migration for a little while.
    Each population undergoes exponential growth after the split.

    Tm is the ratio of split time before migration stops, therefore range 0-1
    """
    nu1, nu2, nu1_fi, nu2_fi, T, m, Tm = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1 * np.exp(np.log(nu1_fi/nu1) * t/T)
    nu2_func = lambda t: nu2 * np.exp(np.log(nu2_fi/nu2) * t/T)
    mu_func = lambda t: m if t < Tm * T else 0

    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=mu_func, m21=mu_func)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def iso_mig_exp(params, ns, pts):
    """
    Split into two populations, with recent migration.
    Each population undergoes exponential growth after the split.

    Tm is the ratio of split time before migration starts, therefore range 0-1
    """
    nu1, nu2, nu1_fi, nu2_fi, T, m, Tm = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: nu1 * np.exp(np.log(nu1_fi/nu1) * t/T)
    nu2_func = lambda t: nu2 * np.exp(np.log(nu2_fi/nu2) * t/T)
    mu_func = lambda t: 0 if t < Tm*T else m

    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=mu_func, m21=mu_func)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def fit_model(data, pts_l, model=no_mig, p0=[10, 10, 1], upper_bound=[100, 100, 10], lower_bound=[1e-2, 1e-2, 0]):
    func = model 
    
    # For no_mig, Parameters are: (nu1, nu2, T)
    # upper_bound = [100, 100, 10]
    # lower_bound = [1e-2, 1e-2, 0]
    
    # This is our initial guess for the parameters, which is somewhat arbitrary.
    # p0 = [10,10,1]
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    
    p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                                  lower_bound=lower_bound)
    
    print('Beginning optimization ************************************************')
    popt = dadi.Inference.opt(p0, data, func_ex, pts_l, 
                              lower_bound=lower_bound,
                              upper_bound=upper_bound,)
                              # verbose=len(p0))
    print('Finshed optimization **************************************************')
    print(popt)
    return popt[0]


def rescale_time(T, theta, num_sites, mu=4.08e-10, gen_per_day=1):
    N_anc = theta / num_sites / mu / 2
    t_gen = T * N_anc
    t_years = t_gen / (365 * gen_per_day)
    return t_years