import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

import config

percid_threshold = 0.1
comp_filter = 200
fig_path = config.ibs_analysis_path / 'figs'
run_summary_path = config.ibs_analysis_path / 'ibs_dat' / f'{config.databatch}__run_year_summary__percid={percid_threshold}__.tsv'

if not fig_path.exists():
    fig_path.mkdir()

run_summary = pd.read_csv(run_summary_path, sep='\t')
run_summary.set_index(['species', 'comp'], inplace=True)
run_summary = run_summary[run_summary['num_comps']>comp_filter]

# rename things so that I don't need to update the plotting code
run_summary['implied_years_log'] = np.log10(run_summary['l99_in_years'])
run_summary['implied_years'] = run_summary['l99_in_years']

allowed_comps = ['Hadza-Tsimane', 'Asia-NorthAmerica', 'Asia-Europe', 'Europe-NorthAmerica', 
                 'Hadza-Hadza', 'Tsimane-Tsimane', 'Asia-Asia', 'Europe-Europe', 'NorthAmerica-NorthAmerica']

####### Violin plot of implied years    #######
fig, ax = plt.subplots(figsize=(5, 4))

sns.violinplot(x='implied_years_log', y='comp', data=run_summary, inner='quart',
               order=allowed_comps, ax=ax, palette='Set2')
comp_median = run_summary.groupby('comp')['implied_years'].median()
ax.set_yticklabels([f'{comp} ({comp_median[comp]:.2f})' for comp in allowed_comps])
# ax.set_xlim(xmin=0)
ax.set_xlabel('L99 years')
ax.set_xticks([3, 4, 5])
ax.set_xticklabels([f'$10^3$', f'$10^4$', f'$10^5$'])
fig.savefig(fig_path / 'violins_log_perc={}.pdf'.format(percid_threshold), bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(figsize=(5, 4))
sns.violinplot(x='implied_years', y='comp', data=run_summary, inner='quart',
               order=allowed_comps, ax=ax, palette='Set2')
comp_median = run_summary.groupby('comp')['implied_years'].median()
ax.set_yticklabels([f'{comp} ({comp_median[comp]:.2f})' for comp in allowed_comps])
ax.set_xlim(xmin=0)
ax.set_xlabel('L99 years')
fig.savefig(fig_path / 'violins_linear_perc={}.pdf'.format(percid_threshold), bbox_inches='tight')


####### Scatter plots #######
run_summary = run_summary.reset_index()
species_summary_df = pd.DataFrame(index=run_summary['species'].unique())

for species, grouped in run_summary.groupby('species'):
    pops = ['Hadza-Hadza', 'Tsimane-Tsimane']
    species_summary_df.loc[species, 'within_HT'] = grouped[grouped['comp'].isin(pops)]['implied_years_log'].mean()
    pops = ['Hadza-Tsimane']
    species_summary_df.loc[species, 'between_HT'] = grouped[grouped['comp'].isin(pops)]['implied_years_log'].mean()
    pops = ['Asia-Asia', 'NorthAmerica-NorthAmerica']
    species_summary_df.loc[species, 'within_CU'] = grouped[grouped['comp'].isin(pops)]['implied_years_log'].mean()
    pops = ['Asia-NorthAmerica']
    species_summary_df.loc[species, 'between_CU'] = grouped[grouped['comp'].isin(pops)]['implied_years_log'].mean()

    pops = ['Europe-Europe', 'NorthAmerica-NorthAmerica']
    species_summary_df.loc[species, 'within_GU'] = grouped[grouped['comp'].isin(pops)]['implied_years_log'].mean()
    pops = ['Europe-NorthAmerica']
    species_summary_df.loc[species, 'between_GU'] = grouped[grouped['comp'].isin(pops)]['implied_years_log'].mean()

    pops = ['Hadza-Nepal']
    species_summary_df.loc[species, 'between_HN'] = grouped[grouped['comp'].isin(pops)]['implied_years_log'].mean()
    pops = ['Nepal-Tsimane']
    species_summary_df.loc[species, 'between_NT'] = grouped[grouped['comp'].isin(pops)]['implied_years_log'].mean()

fig, axes = plt.subplots(1, 3, figsize=(8, 2), dpi=300)
plt.subplots_adjust(wspace=0.5)
s = 5
levels = 4

axes[0].set_aspect('equal', adjustable='box')
axes[1].set_aspect('equal', adjustable='box')
axes[2].set_aspect('equal', adjustable='box')
color = sns.color_palette('pastel', 4)[0]
sns.kdeplot(species_summary_df, ax=axes[0], y='within_HT', x='between_HT', color=color, fill=True, levels=levels)
sns.scatterplot(species_summary_df, ax=axes[0], y='within_HT', x='between_HT', color='k', alpha=0.5, s=s)

color = sns.color_palette('pastel', 4)[1]
sns.kdeplot(species_summary_df, ax=axes[1], y='within_CU', x='between_CU', fill=True, color=color, levels=levels)
sns.scatterplot(species_summary_df, ax=axes[1], y='within_CU', x='between_CU', color='k', alpha=0.5, s=s)

color = sns.color_palette('pastel', 4)[3]
sns.kdeplot(species_summary_df, ax=axes[2], y='within_GU', x='between_GU', fill=True, color=color, levels=levels)
sns.scatterplot(species_summary_df, ax=axes[2], y='within_GU', x='between_GU', color='k', alpha=0.5, s=s)

for ax in axes:
    ax.plot([2, 5], [2, 5], color='k', linestyle='--', linewidth=0.5)
    ax.set_xticks([3, 4, 5])
    ax.set_yticks([3, 4, 5])
    ax.set_xticklabels([f'$10^3$', f'$10^4$', f'$10^5$'])
    ax.set_yticklabels([f'$10^3$', f'$10^4$', f'$10^5$'])
    ax.set_xlabel('Inter-population L99 years')
    ax.set_ylabel('Intra-population L99 years')
    ax.set_xlim(2.5, 5)
    ax.set_ylim(2.5, 5)

axes[0].set_title('Hadza-Tsimane')
axes[1].set_title('Asia-North America')
axes[2].set_title('Europe-North America')
fig.savefig(fig_path / 'scatter_perc={}.pdf'.format(percid_threshold), bbox_inches='tight')
plt.close()
####### global scatter #######
fig, ax = plt.subplots(figsize=(3, 3), dpi=300)
pivot_summary = run_summary.pivot_table(index='species', columns='comp', values='implied_years_log')
sns.kdeplot(pivot_summary,  x='Hadza-Tsimane', y='Asia-NorthAmerica', fill=True, levels=levels)
sns.scatterplot(pivot_summary, x='Hadza-Tsimane', y='Asia-NorthAmerica', color='k', alpha=0.5, s=s)

row = pivot_summary.loc['Escherichia_coli']
plt.scatter(row['Hadza-Tsimane'], row['Asia-NorthAmerica'], color='r', s=10, label='E coli')
plt.legend()

ax.plot([3, 5], [3, 5], color='k', linestyle='--', linewidth=0.5)
ax.set_xticks([3, 4, 5])
ax.set_yticks([3, 4, 5])
ax.set_xticklabels([f'$10^3$', f'$10^4$', f'$10^5$'])
ax.set_yticklabels([f'$10^3$', f'$10^4$', f'$10^5$'])
plt.xlabel('Hadza-Tsimane L99 years')
plt.ylabel('Asia-NorthAmerica L99 years')
plt.savefig(fig_path / 'scatter_global_perc={}.pdf'.format(percid_threshold), bbox_inches='tight')
plt.close()

####### within pop pair plots #######
pivot_summary = run_summary.pivot_table(index='species', columns='comp', values='99_perc_length')
allowed_comps = ['Hadza-Hadza', 'Tsimane-Tsimane', 'Asia-Asia', 'NorthAmerica-NorthAmerica', 'Europe-Europe']
name_map = {
    'Hadza-Hadza': 'Hadza',
    'Tsimane-Tsimane': 'Tsimane',
    'Asia-Asia': 'E. Asia',
    'NorthAmerica-NorthAmerica': 'N. America',
    'Europe-Europe': 'Europe'
}
pivot_summary = pivot_summary[allowed_comps]
pivot_summary.columns = [name_map[col] for col in pivot_summary.columns]

grid = sns.pairplot(data=pivot_summary, corner=True, diag_kind=None, height=1.5)
for i in range(5):
    for j in range(i):
        ax = grid.axes[i, j]
        xs = np.linspace(0, 10000, 100)
        ax.plot(xs, xs, color='k', linestyle='--')
        ax.set_xlim(100, 10000)
        ax.set_ylim(100, 10000)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.setp(ax.get_xticklabels(), fontsize=6)
        plt.setp(ax.get_yticklabels(), fontsize=6)
    grid.axes[i, i].set_visible(False)
grid.savefig(fig_path / 'pairplot_perc={}.pdf'.format(percid_threshold), bbox_inches='tight')
grid.savefig(config.supp_fig_path / 'L99_pairplot.pdf', bbox_inches='tight')

plt.close()
