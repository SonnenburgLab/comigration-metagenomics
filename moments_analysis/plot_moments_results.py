import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import config

pops = ['Hadza', 'Tsimane']
model_name = 'split_mig'
# pops = ['MetaHIT', 'HMP']
popstr = ''.join(pops)
moments_results = pd.read_csv(config.project_path / 'moments_analysis' / f'moments_results_{popstr}_{model_name}.csv', index_col=0)
moments_figure_path = config.fig_path / '241118_moments'
moments_figure_path.mkdir(exist_ok=True)

print("Print {} passed genome filter".format(moments_results.shape[0]))

# filter to species with clade structure
passed_species = moments_results[~moments_results['if_clade']]
print("Print {} passed clade filter".format(passed_species.shape[0]))

# filter to species with low uncertainty
uncert_ratio = passed_species['Tsplit_uncert'] / passed_species['Tsplit']

# if the uncertainty on the split time is too large, drop the species
passed = (uncert_ratio < 0.5) & (np.isnan(uncert_ratio)==False)
passed_species = passed_species[passed]
print("Print {} passed uncertainty filter".format(sum(passed)))

passed_species_after_uncert = passed_species.copy()
passed_species = passed_species[passed_species['mean_abs_resid'] < 1.5]
print("Print {} passed residual filter".format(passed_species.shape[0]))


################################################################
# first plot residual distribution
################################################################
fig, axes = plt.subplots(1, 2, figsize=(5, 3),)  # Adjust width for better visibility
plt.subplots_adjust(wspace=0.5)

# Panel 1: Mean residual size
sns.boxplot(
    data=moments_results, 
    y='mean_abs_resid', 
    x='if_clade', 
    hue='if_clade', 
    ax=axes[0]
)
num_false = moments_results[moments_results['if_clade'] == False].shape[0]
num_true = moments_results[moments_results['if_clade'] == True].shape[0]
axes[0].set_xticks([0, 1])
axes[0].set_xticklabels([f'False\n(n={num_false})', f'True\n(n={num_true})'])
axes[0].legend().remove()
axes[0].set_xlabel('')
# axes[0].set_title('Mean Residual Size', fontsize=10)
axes[0].set_ylabel('Mean residual size')

# Panel 2: Split time (log scale)
sns.boxplot(
    data=moments_results, 
    y='Tsplit', 
    x='if_clade', 
    hue='if_clade', 
    ax=axes[1]
)
axes[1].set_yscale('log')  # Log scale for the y-axis
axes[1].set_xticks([0, 1])
axes[1].set_xticklabels([f'False\n(n={num_false})', f'True\n(n={num_true})'])
axes[1].legend().remove()
axes[1].set_xlabel('')
# axes[1].set_title('Split Time (Years)', fontsize=10)
axes[1].set_ylabel('Split time (years)')

axes[0].set_xlabel("If clade structure")
axes[1].set_xlabel("If clade structure")
plt.savefig(moments_figure_path / f'residual_size_by_clades_{popstr}_{model_name}.pdf', bbox_inches='tight')
plt.close()

################################################################
# second plot residual among species with no clades
################################################################
plt.figure(figsize=(6, 2))

to_plot = passed_species_after_uncert
to_plot['passed'] = to_plot['mean_abs_resid'] < 1.5

sns.barplot(data=to_plot, x='species_names', y='mean_abs_resid', hue='passed',
            linewidth=1.5, edgecolor=".5", palette=['white', 'tab:blue'])
plt.legend().remove()
plt.axhline(1.5, color='k', linestyle='--')
_ = plt.xticks(rotation=90)
plt.ylabel('Mean residual size')
plt.xlabel('')
plt.savefig(moments_figure_path / f'residual_size_by_species_{popstr}_{model_name}.pdf', bbox_inches='tight')
plt.close()


################################################################
# third plot split time values
################################################################
plt.figure(figsize=(6, 2))

# load human parameters
human_model = pd.read_csv(config.intermediate_data_path / 'MedinaMunoz_model.csv', index_col=0)
human_model.iloc[:, 1:] = human_model.iloc[:, 1:].astype(float)

plot_results = passed_species

# first plot bacterial results
bacterial_locs = np.array(range(plot_results.shape[0])) + 3
plt.errorbar(bacterial_locs, plot_results['Tsplit'] / 1000, yerr=2*plot_results['Tsplit_uncert'] / 1000, fmt='o', color='tab:blue')
bacterial_names = plot_results.index.to_list()

# next plot human results
loc = 0
values = ['TB', 'TN']
colors = ['tab:blue', 'tab:blue']
human_locs = []
for pair in zip(values, colors):
    value = pair[0]
    color = pair[1]
    # plt.errorbar(loc-0.2, human_model.loc[value, 'intronic_value'], yerr=2*human_model.loc[value, 'intronic_SE'], fmt='o', color=color)
    # plt.errorbar(loc+0.2, human_model.loc[value, 'intergenic_value'], yerr=2*human_model.loc[value, 'intergenic_SE'], fmt='o', color=color)
    human_locs.append(loc)
    # plotting only inferred values using synonymous sites
    plt.errorbar(loc, human_model.loc[value, 'syn_value'], yerr=2*human_model.loc[value, 'syn_SE'], fmt='o', color=color)
    plt.axhline(human_model.loc[value, 'syn_value'], color='grey', linestyle='--')
    loc += 1

plt.xticks(human_locs + list(bacterial_locs), ['Out-of-Africa split', 'Asia-Americas split'] + bacterial_names, 
           rotation=-45, ha='left')
plt.axvspan(-1, (max(human_locs) + min(bacterial_locs))/2, color='gray', alpha=0.2)

plt.ylabel('Time (years ago)')
plt.yscale('log')
plt.ylim(10, 1e3)
plt.yticks([10, 100, 1000])
plt.gca().set_yticklabels(['10,000', '100,000', '1,000,000'])
plt.xlim(-1, max(bacterial_locs)+1)
plt.savefig(moments_figure_path / f'split_time_values_{popstr}_{model_name}.pdf', bbox_inches='tight')
plt.close()

# finally, save the moments results to supp table
# renaming some columns
plot_results = plot_results.rename(columns={'uncert_nu1': 'nu1_std', 'uncert_nu2': 'nu2_std', 'uncert_T': 'T_std', 'uncert_m': 'm_std','uncert_theta': 'theta_std', 'num_Hadza': 'num_Hadza_MAGs', 'num_Tsimane': 'num_Tsimane_MAGs', 'proj_Hadza': 'Hadza_projection', 'proj_Tsimane': 'Tsimane_projection', 'syn_genome_length': 'num_total_syn_sites', 'Tsplit': 'Tsplit (yr)', 'Tsplit_uncert': 'Tsplit_std (yr)'})

columns_to_keep = ['nu1', 'nu1_std', 'nu2', 'nu2_std', 'T', 'T_std', 'm', 'm_std', 'theta', 'theta_std', 'num_Hadza_MAGs', 'num_Tsimane_MAGs', 'Hadza_projection', 'Tsimane_projection', 'num_total_syn_sites', 'num_sites_passing_proj', 'Tsplit (yr)', 'Tsplit_std (yr)']
plot_results = plot_results[columns_to_keep]
plot_results.to_csv(config.supp_table_path / f'supp_moments_results_{popstr}_{model_name}.tsv', sep='\t')