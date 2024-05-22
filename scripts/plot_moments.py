import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

import moments
from utils import moments_utils, metadata_utils

data_batch = '240509'
focal_pops = ['Hadza', 'Tsimane']

model_name = 'split_no_mig'
param_names = ['nu1', 'nu2', 'T']
model_func = moments_utils.split_no_mig

# model_name = 'split_mig'
# model_func = moments.Demographics2D.split_mig
# param_names = ['nu1', 'nu2', 'T', 'm']

output_pdf = f'figs/moments_sfs__{data_batch}__{model_name}.pdf'

moments_results = pd.read_csv(f'moments_out/{data_batch}__{model_name}.csv', index_col=False)
moments_results.set_index('species', inplace=True)

metahelper = metadata_utils.MetadataHelper(data_batch=data_batch)

with PdfPages(output_pdf) as pdf:
    for species, row in moments_results.iterrows():
        logging.info(f"Processing {species}")

        data, _ = res = moments_utils.load_SFS(species, metahelper, focal_pops=focal_pops)
        data.mask[1, :] = True
        data.mask[:, 1] = True

        opt_params = row[param_names]

        plt.figure(figsize=(8,6))
        moments.Plotting.plot_2d_comp_multinom(
            model_func(opt_params, data.sample_sizes), data, resid_range=5)

        plt.gcf().text(0.5, 1.1, species, ha='center', va='top')
        param_string = ', '.join([f"{param}: {val:.2f}" for param, val in zip(param_names, opt_params)])
        plt.gcf().text(0.5, 1., param_string, ha='center', va='top')

        pdf.savefig(bbox_inches='tight')
        plt.close()