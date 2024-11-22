import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import config

# pops = ['Hadza', 'Tsimane']
model_name = 'split_no_mig'
pops = ['MetaHIT', 'HMP']
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
plt.figure(figsize=(1.5, 3))
sns.boxplot(moments_results, y='mean_abs_resid', x='if_clade', hue='if_clade')
num_false = sum(moments_results['if_clade']==False)
num_true = sum(moments_results['if_clade']==True)
plt.xticks([0, 1], ['False\n(n={})'.format(num_false), 'True\n(n={})'.format(num_true)])
plt.legend().remove()
plt.xlabel('')
plt.title('Potential clade structure?', fontsize=10)
plt.ylabel('Mean residual size')
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