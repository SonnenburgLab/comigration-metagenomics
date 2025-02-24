import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import moments

from utils import moments_utils
import config

moment_res_path = config.supp_table_path / 'supp_moments_results_HadzaTsimane_split_mig.tsv'
moments_results = pd.read_csv(moment_res_path, sep='\t', index_col=0)

font_size = 3
mpl.rcParams['axes.titlesize'] = font_size
mpl.rcParams['axes.labelsize'] = font_size
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size
mpl.rcParams['legend.fontsize'] = font_size

def plot_residual(ax, data, model, mask_singletons=False):
    if mask_singletons:
        data.mask[1, :] = True
        data.mask[:, 1] = True
    vmin = 1
    resid_range = 5
    pops = data.pop_ids
    resid = moments.Inference.linear_Poisson_residual(model, data, mask=vmin)
    moments.Plotting.plot_2d_resid(
        resid,
        resid_range,
        pop_ids=pops,
        ax=ax,
    )
    return resid

passed_species = moments_results.index.to_list()

####################################################################################################
# Plot SFS and residuals for the first 8 species
####################################################################################################
species_to_plot = passed_species[:8]
fig, axes = plt.subplots(len(species_to_plot), 6, figsize=(6.5, 8))
plt.subplots_adjust(hspace=0.7, wspace=0.5)

for i, species in enumerate(species_to_plot):
    proj = moments_results.loc[species, ['Hadza_projection', 'Tsimane_projection']].astype(int).to_list()
    data = moments_utils.load_SFS_projection(species, focal_pops=['Hadza', 'Tsimane'])
    model = moments_utils.prep_model(moments_results.loc[species], moments.Demographics2D.split_mig, data,
                                     model_name='split_mig')
    f_h = data.marginalize([1])
    f_t = data.marginalize([0])
    model_h = model.marginalize([1])
    model_t = model.marginalize([0])

    axes[i, 0].plot(f_t, 'o', color='tab:blue', label='Data', alpha=0.8, markersize=3, markerfacecolor='none')
    axes[i, 0].plot(model_t, '-o', color='tab:orange', label='Model', alpha=0.8, markersize=1.5, linewidth=0.5)
    axes[i, 0].set_title(f'Tsimane', pad=2)
    axes[i, 0].set_yscale('log')
    axes[i, 0].set_xlabel('Minor allele count', labelpad=0.5)
    axes[i, 0].set_ylabel(species + '\n\nNum observed', labelpad=0.5)

    axes[i, 1].plot(f_h, 'o', color='tab:blue', label='Data', alpha=0.8, markersize=3, markerfacecolor='none')
    axes[i, 1].plot(model_h, '-o', color='tab:orange', label='Model', alpha=0.8, markersize=1.5, linewidth=0.5)
    axes[i, 1].set_xlabel('Minor allele count', labelpad=0.5)
    axes[i, 1].set_title(f'Hadza', pad=2)
    axes[i, 1].set_yscale('log')


    axes[i, 0].axvspan(0, 1.5, color='gray', alpha=0.5)
    axes[i, 0].set_xlim(xmin=0)
    axes[i, 1].axvspan(0, 1.5, color='gray', alpha=0.5)
    axes[i, 1].set_xlim(xmin=0)

    moments.Plotting.plot_single_2d_sfs(data, vmin=1, ax=axes[i, 2], cmap='viridis')
    axes[i, 2].set_xlabel('Tsimane', labelpad=0)
    axes[i, 2].set_ylabel('Hadza', labelpad=0)
    axes[i, 2].set_title('Data', pad=2)
    moments.Plotting.plot_single_2d_sfs(model, vmin=1, ax=axes[i, 3], cmap='viridis',
                                        pop_ids=['Hadza', 'Tsimane'])
    axes[i, 3].set_xlabel('Tsimane', labelpad=0)
    axes[i, 3].set_ylabel('Hadza', labelpad=0)
    axes[i, 3].set_title('Model', pad=2)

    resid = plot_residual(axes[i, 4], data, model, mask_singletons=True)
    axes[i, 4].set_xlabel('Tsimane', labelpad=0)
    axes[i, 4].set_ylabel('Hadza', labelpad=0)
    axes[i, 4].set_title(r'$\mathrm{Resid}=(\mathrm{Model}-\mathrm{Data})/\sqrt{\mathrm{Model}}$', pad=2)

    flatresid = np.compress(np.logical_not(resid.mask.ravel()), resid.ravel())
    bins = np.linspace(-5, 5, 20)
    axes[i, 5].hist(flatresid, bins=bins, color='tab:blue', density=True)
    axes[i, 5].set_xlim(-5, 5)
    axes[i, 5].set_yticks([])
    axes[i, 5].set_title('Resid histogram', pad=2)
axes[0, 0].legend()

fig.savefig(config.supp_fig_path / 'moments_sfs_1.pdf', dpi=300, bbox_inches='tight')
plt.close()

####################################################################################################
# Plot SFS and residuals for the remaining species
####################################################################################################
species_to_plot = passed_species[8:]
fig, axes = plt.subplots(len(species_to_plot), 6, figsize=(6.5, len(species_to_plot)))
plt.subplots_adjust(hspace=0.7, wspace=0.5)

for i, species in enumerate(species_to_plot):
    proj = moments_results.loc[species, ['Hadza_projection', 'Tsimane_projection']].astype(int).to_list()
    data = moments_utils.load_SFS_projection(species, focal_pops=['Hadza', 'Tsimane'])
    model = moments_utils.prep_model(moments_results.loc[species], moments.Demographics2D.split_mig, data,
                                     model_name='split_mig')
    f_h = data.marginalize([1])
    f_t = data.marginalize([0])
    model_h = model.marginalize([1])
    model_t = model.marginalize([0])

    axes[i, 0].plot(f_t, 'o', color='tab:blue', label='Data', alpha=0.8, markersize=3, markerfacecolor='none')
    axes[i, 0].plot(model_t, '-o', color='tab:orange', label='Model', alpha=0.8, markersize=1.5, linewidth=0.5)
    axes[i, 0].set_title(f'Tsimane', pad=2)
    axes[i, 0].set_yscale('log')
    axes[i, 0].set_xlabel('Minor allele count', labelpad=0.5)
    axes[i, 0].set_ylabel(species + '\n\nNum observed', labelpad=0.5)

    axes[i, 1].plot(f_h, 'o', color='tab:blue', label='Data', alpha=0.8, markersize=3, markerfacecolor='none')
    axes[i, 1].plot(model_h, '-o', color='tab:orange', label='Model', alpha=0.8, markersize=1.5, linewidth=0.5)
    axes[i, 1].set_xlabel('Minor allele count', labelpad=0.5)
    axes[i, 1].set_title(f'Hadza', pad=2)
    axes[i, 1].set_yscale('log')


    axes[i, 0].axvspan(0, 1.5, color='gray', alpha=0.5)
    axes[i, 0].set_xlim(xmin=0)
    axes[i, 1].axvspan(0, 1.5, color='gray', alpha=0.5)
    axes[i, 1].set_xlim(xmin=0)

    moments.Plotting.plot_single_2d_sfs(data, vmin=1, ax=axes[i, 2], cmap='viridis')
    axes[i, 2].set_xlabel('Tsimane', labelpad=0)
    axes[i, 2].set_ylabel('Hadza', labelpad=0)
    axes[i, 2].set_title('Data', pad=2)
    moments.Plotting.plot_single_2d_sfs(model, vmin=1, ax=axes[i, 3], cmap='viridis',
                                        pop_ids=['Hadza', 'Tsimane'])
    axes[i, 3].set_xlabel('Tsimane', labelpad=0)
    axes[i, 3].set_ylabel('Hadza', labelpad=0)
    axes[i, 3].set_title('Model', pad=2)

    resid = plot_residual(axes[i, 4], data, model, mask_singletons=True)
    axes[i, 4].set_xlabel('Tsimane', labelpad=0)
    axes[i, 4].set_ylabel('Hadza', labelpad=0)
    axes[i, 4].set_title(r'$\mathrm{Resid}=(\mathrm{Model}-\mathrm{Data})/\sqrt{\mathrm{Model}}$', pad=2)

    flatresid = np.compress(np.logical_not(resid.mask.ravel()), resid.ravel())
    bins = np.linspace(-5, 5, 20)
    axes[i, 5].hist(flatresid, bins=bins, color='tab:blue', density=True)
    axes[i, 5].set_xlim(-5, 5)
    axes[i, 5].set_yticks([])
    axes[i, 5].set_title('Resid histogram', pad=2)
axes[0, 0].legend()

fig.savefig(config.supp_fig_path / 'moments_sfs_2.pdf', dpi=300, bbox_inches='tight')
plt.close()