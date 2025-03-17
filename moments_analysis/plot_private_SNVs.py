import pandas as pd
import moments
import matplotlib.pyplot as plt

from utils import moments_utils
import config

def compute_and_save_likelihoods(sfs_batch, passed_species):
    save_path = config.moments_path / 'moments_dat' / f'{sfs_batch}_private_SNV_proba'
    if not save_path.exists():
        save_path.mkdir(parents=True)
    for species in passed_species.index:
        print(species)
        sfs_dat = moments_utils.load_sfs_counts(species, sfs_folder=config.sfs_path / sfs_batch)
        sfs_dat = moments_utils.compute_well_mixed_likelihoods(sfs_dat)

        sample_sizes = sfs_dat['Hadza_tot'].max(), sfs_dat['Tsimane_tot'].max()

        # prepare inferred demographic model
        dem_model = moments_utils.prep_model_unscaled(passed_species.loc[species], moments.Demographics2D.split_mig, sample_sizes=sample_sizes)
        # prepare high migration model
        high_mig = passed_species.loc[species].copy()
        high_mig['m'] = 100
        model_high_mig = moments_utils.prep_model_unscaled(high_mig, moments.Demographics2D.split_mig, sample_sizes=sample_sizes)
        model_high_mig = model_high_mig / model_high_mig.sum()

        model_probs = moments_utils.compute_exclusive_likelihood_model(sfs_dat, dem_model)
        sfs_dat = sfs_dat.join(model_probs, rsuffix='_inferred')

        model_probs = moments_utils.compute_exclusive_likelihood_model(sfs_dat, model_high_mig)
        sfs_dat = sfs_dat.join(model_probs, rsuffix='_high_mig')
        sfs_dat.to_csv(save_path / f'{species}_exclusive_snv_likelihoods.csv')


if __name__ == '__main__':
    sfs_batch = config.sfs_batch
    model_name = 'split_mig'
    pops = ['Hadza', 'Tsimane']
    # moment_results = pd.read_csv(config.moments_path / 'moments_dat' / f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}_full.csv')
    # passed_species = pd.read_csv(config.intermediate_data_path / 'fig4_data.csv', index_col=0)
    passed_species = pd.read_csv(config.moments_path / 'moments_dat' / f'{sfs_batch}__{model_name}__{pops[0]}__{pops[1]}_cleaned.csv', index_col=0)

    # only need to run once
    # compute_and_save_likelihoods(sfs_batch, passed_species)

    species_results = {}
    for species in passed_species.index:
        sfs_dat = pd.read_csv(config.moments_path / 'moments_dat' / f'{sfs_batch}_private_SNV_proba' / f'{species}_exclusive_snv_likelihoods.csv', index_col=0)
        species_results[species] = sfs_dat

    ################################################################
    # plot the "supplementary panels"
    ################################################################
    # color_null = '#008080'
    # color_inferred = '#FF6F61'
    # color_high_mig = '#66B2B2'
    # color_data = '#2F4F4F'

    # match Matt's color choices
    color_data = '#f1be11'
    color_null = '#10193b'
    color_inferred = '#f9ec37'

    fig, axes = plt.subplots(len(passed_species)//4 + 1, 4, figsize=(9, 7))
    plt.subplots_adjust(hspace=0.5, wspace=0.4)
    axes_flat = axes.flatten()

    for i, species in enumerate(passed_species.index):
        sfs_dat = species_results[species]
        obs_prob = sfs_dat.groupby('Total_alt')['Is_exclusive'].mean()
        # well-mixed model
        null_model_prob = sfs_dat.groupby('Total_alt')['Exclusive_snv_likelihood'].mean()
        # inferred model
        model_prob = sfs_dat.groupby('Total_alt')['Model_exclusive_likelihood'].mean()
        # high migration model
        model_high_mig_prob = sfs_dat.groupby('Total_alt')['Model_exclusive_likelihood_high_mig'].mean()

        ax = axes_flat[i]
        ax.plot(obs_prob, label='Observed', marker='o', linestyle='', color=color_data, alpha=0.5, zorder=4)
        ax.plot(null_model_prob, linewidth=2, label='No Isolation', color=color_null)
        # ax.plot(model_high_mig_prob, label='High Migration\n$N_em=100$', linestyle='--', color=color_high_mig)
        ax.plot(model_prob, linewidth=2, label='Inferred model', color=color_inferred)

        ax.set_xlim(1, 25)
        ax.set_ylim(1e-3, 1)
        ax.set_yscale('log')

        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_title(species, fontsize=8)

    # for i in range(3):
    #     for j in range(1, 3):
    #         axes[i, j].set_ylabel('')
    # for i in range(2):
    #     for j in range(3):
    #         axes[i, j].set_xlabel('')
    axes[-1, -1].set_visible(False)
    axes[-1, -2].legend(bbox_to_anchor=(1.2, 0.5), loc='center left', frameon=False)

    # Add shared X and Y labels
    fig.text(0.5, 0.04, 'Total alt allele counts', ha='center', va='center', fontsize=12)
    fig.text(0.04, 0.5, 'Fraction of population-specific SNVs', ha='center', va='center', rotation='vertical', fontsize=12)

    # plt.savefig(config.fig_path / '250129_exclusive_snv_likelihoods.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(config.supp_fig_path / 'private_SNV_proba.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    ################################################################
    # plot the "main panel"
    # No longer needed -- Matt's plotting this in R
    ################################################################
    # example_species_name = 'UBA738_sp003522945'
    # sfs_dat = species_results[example_species_name]
    # avg_tot = sfs_dat['Hadza_tot'].mean() + sfs_dat['Tsimane_tot'].mean()

    # fig, ax = plt.subplots(figsize=(2.5, 2))

    # # overall logic, for each alt allele count, compute the fraction of SNVs that are exclusive to one population for observed data
    # # For simulated models, compute the mean likelihood of exclusive SNVs for each alt allele count

    # # well-mixed model; computed from binomial sampling likelihood
    # null_model_prob = sfs_dat.groupby('Total_alt')['Exclusive_snv_likelihood'].mean()
    # xs = null_model_prob.index.values / avg_tot

    # # Is_exclusive will be 0/1 for each SNV, so the mean is the fraction of SNVs that are exclusive to one population
    # obs_prob = sfs_dat.groupby('Total_alt')['Is_exclusive'].mean()
    # xs = obs_prob.index.values / avg_tot
    # ax.plot(xs, obs_prob, label='Observed', marker='o', linestyle='', color=color_data, alpha=0.3, zorder=10)

    # colors = [color_null, color_inferred]
    # true_split = passed_species.loc[example_species_name, 'Tsplit (yr)']
    # labels = ['No split', 'T={}yr\n(Inferred model)'.format(int(true_split))]
    # col_vals = ['Exclusive_snv_likelihood', 'Model_exclusive_likelihood']
    # for i in range(2):
    #     model_prob = sfs_dat.groupby('Total_alt')[col_vals[i]].mean()
    #     xs = model_prob.index.values / avg_tot
    #     ax.plot(xs, model_prob, label=labels[i], color=colors[i])

    # # ax.set_xlim(1, 20)
    # ax.set_xlim(0, 0.2)
    # ax.set_ylim(1e-3, 1)
    # ax.set_yscale('log')
    # ax.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=False)

    # ax.set_xlabel("SNV frequency")
    # ax.set_ylabel("Fraction of SNVs \nin only one pop")
    # ax.set_title(example_species_name)
    # plt.savefig(config.fig_path / '20250313_private_snv_proba_inferred_example.pdf', dpi=300, bbox_inches='tight')