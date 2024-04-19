import dadi
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from matplotlib.backends.backend_pdf import PdfPages
from utils import dadi_utils

output_pdf = 'figs/dadi_output_20240319_masked_100-.pdf'
output_csv = 'figs/dadi_output_20240319_masked100-.csv'
snp_folder = 'dat/240318_dadi_snps/'
species_df = pd.read_csv(snp_folder + '_metadata.csv')
mask_singletons = True

# find the top species with a lot of samples
species_df['Num Total'] = species_df['Num Hadza'] + species_df['Num Tsimane']
species_df.set_index('Species', inplace=True)
species_df = species_df[species_df['Num Total']<100]
print("Total number of species:", species_df.shape[0])

# fit species
with PdfPages(output_pdf) as pdf:
    for species, row in species_df.iterrows():
        num_hadza = row['Num Hadza']
        num_ts = row['Num Tsimane']
        print(species)
        if num_hadza < 5 or num_ts < 5:
            print("Skipping species with less than 5 samples in either population")
            continue
        
        snp_file = os.path.join(snp_folder, species + '.snps.txt')
        dd = dadi.Misc.make_data_dict(snp_file)
        
        pop_ids=["Hadza", "Tsimane"]
        # TODO: this is a bit arbitrary, but keep sites that are covered in >85% of samples
        # if lower than 5, I really don't think we can get anything out of the data
        h_proj = max(5, int(num_hadza * 0.85))
        ts_proj = max(5, int(num_ts * 0.85))
        proj = [h_proj, ts_proj]
    
        data = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)
        if mask_singletons:
            # masking the singletons
            data.mask[1, :] = True
            data.mask[:, 1] = True
        ns = data.sample_sizes
        
        min_pts = int(max(num_hadza, num_ts))  # arbitrary choice to ensure grid is finer than samples
        pts_l = [min_pts, int(min_pts * 1.5), int(min_pts * 2)]

        params = dadi_utils.fit_model(data=data, model=dadi_utils.no_mig, pts_l=pts_l)
        
        func_ex = dadi.Numerics.make_extrap_log_func(dadi_utils.no_mig)

        model = func_ex(params, ns, pts_l)
        # Likelihood of the data given the model AFS.
        ll_model = dadi.Inference.ll_multinom(model, data)
        # The optimal value of theta given the model.
        theta = dadi.Inference.optimal_sfs_scaling(model, data)
        print('Optimal value of theta: {0}'.format(theta))
    
        plt.figure(figsize=(8,6))
        dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=5,
                                            pop_ids =('Hazda','Tsimane'), show=False)

        species_df.loc[species, ['nu1', 'nu2', 'T split']] = params
        species_df.loc[species, 'theta'] = theta
        
        plt.gcf().text(0.5, 1.1, species, ha='center', va='top')
        plt.gcf().text(0.5, 1., "nu1={0:.2f}; nu2={1:.2f}; T={2:.3f}".format(params[0], params[1], params[2]), ha='center', va='top')
        
        pdf.savefig(bbox_inches='tight')
        plt.close()
        print()

species_df.to_csv(output_csv, index=True)