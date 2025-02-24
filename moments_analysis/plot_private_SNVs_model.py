import moments
import pandas as pd
import matplotlib.pyplot as plt

from utils import moments_utils
import config

def compute_model():
    # choose the species observed counts
    species = 'UBA738_sp003522945'

    # put in model parameters by hand
    test_sfs_dat = moments_utils.load_sfs_counts(species)
    test_sfs_dat = moments_utils.compute_well_mixed_likelihoods(test_sfs_dat)
    sample_sizes = test_sfs_dat['Hadza_tot'].max(), test_sfs_dat['Tsimane_tot'].max()

    test_ms = [0, 1, 5]
    test_Ts = [0, 0.1, 0.2, 0.5, 1, 2]
    for i, m in enumerate(test_ms):
        test_model_params = pd.DataFrame({'nu1': [2], 'nu2': [2], 'T': [1], 'm': [m], 'theta': [1000]}, index=['Test'], dtype=float).loc['Test']

        for j, T in enumerate(test_Ts):
            # prepare inferred demographic model
            test_model_params['T'] = T
            test_model = moments_utils.prep_model_unscaled(test_model_params, moments.Demographics2D.split_mig, sample_sizes=sample_sizes)
            test_model = test_model / test_model.sum()

            model_probs = moments_utils.compute_exclusive_likelihood_model(test_sfs_dat, test_model)
            test_sfs_dat = test_sfs_dat.join(model_probs, rsuffix=f'_{i}_{j}')

    sfs_dat = test_sfs_dat
    sfs_dat.to_csv(config.intermediate_data_path / '250129_private_SNVs_model.csv')

if __name__ == '__main__':
    # compute_model()
    sfs_dat = pd.read_csv(config.intermediate_data_path / '250129_private_SNVs_model.csv', index_col=0)

    test_ms = [0, 1, 5]
    test_Ts = [0, 0.1, 0.2, 0.5, 1, 2]
    avg_tot = sfs_dat['Hadza_tot'].mean() + sfs_dat['Tsimane_tot'].mean()
    fig, axes = plt.subplots(1, 3, figsize=(8, 2))
    plt.subplots_adjust(wspace=0.3)

    for i, m in enumerate(test_ms):
        ax = axes[i]
        for j in range(len(test_Ts))[::-1]:
            if j == 0:
                label = 'Model_exclusive_likelihood'
            else:
                label = f'Model_exclusive_likelihood_{i}_{j}'
            # inferred model
            model_prob = sfs_dat.groupby('Total_alt')[label].mean()
            xs = model_prob.index.values / avg_tot
            ax.plot(xs, model_prob, label='T={}'.format(test_Ts[j]), color='C{}'.format(j))

        # ax.set_xlim(1, 20)
        ax.set_xlim(0, 0.2)
        ax.set_ylim(1e-3, 1)
        ax.set_yscale('log')

        ax.set_xlabel("SNV frequency")
        ax.set_title('m={}'.format(m))
    axes[-1].legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
    axes[0].set_ylabel("Fraction of SNVs\nin only one pop")

    # plt.savefig(config.fig_path / '250129_private_SNVs_model.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(config.supp_fig_path/ 'private_SNV_proba_model.pdf', dpi=300, bbox_inches='tight')

