"""
Utilities for manipulating SNV catalogs.
"""

import os
import pandas as pd
import numpy as np
import config

def load_degeneracy(degeneracy_file):
    degen_sites = pd.read_csv(degeneracy_file, delimiter='\t')
    # degen_sites.columns = ['contig', 'location', 'degeneracy', 'ref']
    degen_sites.columns = ['Contig', 'Site', 'Degeneracy', 'Base']
    degen_sites.set_index(['Contig', 'Site'], inplace=True)
    return degen_sites


def load_snvs(snv_file):
    snv_catalog = pd.read_csv(snv_file, delimiter='\t')
    snv_catalog['Contig'] = snv_catalog['snp_id'].apply(lambda x: x.split('||')[0])
    snv_catalog['Site'] = snv_catalog['snp_id'].apply(lambda x: int(x.split('||')[1]))
    snv_catalog['ref'] = snv_catalog['snp_id'].apply(lambda x: x.split('||')[-2])
    snv_catalog['alt'] = snv_catalog['snp_id'].apply(lambda x: x.split('||')[-1])
    snv_catalog.set_index(['Contig', 'Site'], inplace=True)
    return snv_catalog


def load_core_gene_mask(core_gene_file):
    # check if file is empty
    try:
        core_gene_df = pd.read_csv(core_gene_file, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        return None
    core_gene_df.columns = ['Contigs', 'Site']
    core_gene_index = core_gene_df.set_index(['Contigs', 'Site']).index
    return core_gene_index

def load_dadi_snps_metadata(batch_name):
    metadata_path = os.path.join(config.dadi_dat_path, '{}__metadata.csv'.format(batch_name))
    return pd.read_csv(metadata_path, index_col=0)

def get_population_accession_identifier(pop_name):
    if pop_name == 'Hadza':
        return 'Hadza'
    elif pop_name == 'Tsimane':
        return 'TS'
    elif pop_name == 'HMP':
        return 'SRS'
    elif pop_name == 'MetaHIT':
        return 'ERS'
    else:
        raise ValueError('Population not recognized')
    
def sample_name_to_population(sample_name):
    if 'Hadza' in sample_name:
        return 'Hadza'
    elif 'TS' in sample_name:
        return 'Tsimane'
    elif 'SRS' in sample_name:
        return 'HMP'
    elif 'ERS' in sample_name:
        return 'MetaHIT'
    else:
        raise ValueError('Population not recognized')

def sample_number_from_header(snv_file):
    """
    Extract the number of samples from the header of a snv file.
    """
    # read first line to get the column names
    with open(snv_file) as f:
        first_line = f.readline().strip()
        samples = first_line.split('\t')[1:]
        num_hadza, num_ts = sample_names_to_H_T(samples)
    return num_hadza, num_ts

def sample_names_to_H_T(samples):
    num_hadza = sum(['Hadza' in col for col in samples])
    num_ts = sum(['TS' in col for col in samples])
    return num_hadza, num_ts

def sample_names_to_HMP(samples):
    # TODO: eventually use the true accession table
    num_hmp = sum(['SRS' in col for col in samples])
    return num_hmp

def sample_names_to_metahit(samples):
    # TODO: eventually use the true accession table
    num_meta = sum(['ERS' in col for col in samples])
    return num_meta

def load_and_filter_snv_catalog(snv_file, degeneracy_file, core_gene_file):
    """
    Returns only 1D and 4D sites in the core genome.
    Shape: (n_snps, n_samples)
    """

    degen_sites = load_degeneracy(degeneracy_file)
    snv_catalog = load_snvs(snv_file)
    core_gene_index = load_core_gene_mask(core_gene_file)

    degen_1D = degen_sites[degen_sites['Degeneracy'] == 1]
    degen_4D = degen_sites[degen_sites['Degeneracy'] == 4]

    core_1D = snv_catalog.index.isin(degen_1D.index) & snv_catalog.index.isin(core_gene_index)
    core_4D = snv_catalog.index.isin(degen_4D.index) & snv_catalog.index.isin(core_gene_index)

    syn_snvs = snv_catalog.loc[core_4D, :].copy()
    syn_snvs.set_index(['ref', 'alt'], append=True, inplace=True)
    syn_snvs.drop(['snp_id'], axis=1, inplace=True)
    nonsyn_snvs = snv_catalog.loc[core_1D, :].copy()
    nonsyn_snvs.set_index(['ref', 'alt'], append=True, inplace=True)
    nonsyn_snvs.drop(['snp_id'], axis=1, inplace=True)

    # saving previous approach for now, copying degeneracy and filter;
    # also took care of duplicate sites in degenracy file

    # # only coding sites are in the degeneracy file
    # coding_mask = snv_catalog.index.isin(degen_sites.index)
    # # remove duplicate sites
    # # these are likely caused by overlapping genes, so two annotations
    # dup_sites = degen_sites.index[degen_sites.index.duplicated()]
    # dup_mask = snv_catalog.index.isin(dup_sites)
    # filtered_snps = snv_catalog.loc[coding_mask & (~dup_mask), :].copy()
    # # copy degeneracy info to dataframe
    # filtered_snps['Degeneracy'] = degen_sites[~degen_sites.index.duplicated()]['Degeneracy']
    # # put all necessary columns in the index
    # filtered_snps.set_index(['ref', 'alt'], append=True, inplace=True)
    # # for most use cases, only 4D and 1D sites are useful
    # syn_snps = filtered_snps.loc[filtered_snps['Degeneracy']==4, :].copy()
    # nonsyn_snps = filtered_snps.loc[filtered_snps['Degeneracy']==1, :].copy()
    # syn_snps.drop(['snp_id', 'Degeneracy'], axis=1, inplace=True)
    # nonsyn_snps.drop(['snp_id', 'Degeneracy'], axis=1, inplace=True)
    return syn_snvs, nonsyn_snvs


def get_core_genome_length(degeneracy_file, core_gene_file):
    """
    Return both the total length of the core genome and the length of the 4D core genome.
    """
    degen_sites = load_degeneracy(degeneracy_file)
    core_gene_index = load_core_gene_mask(core_gene_file)
    if core_gene_index is None:
        return 0, 0
    core_4D = degen_sites.index.isin(core_gene_index)
    return core_gene_index.shape[0], core_4D.sum()

def split_snvs(snv_catalog, pops=['Hadza', 'TS']):
    """
    Split the snv_catalog into multiple dataframes, one for each population.
    """
    pop_snvs = []
    for pop in pops:
        pop_mask = np.array([pop in x for x in snv_catalog.columns])
        pop_snvs.append(snv_catalog.loc[:, pop_mask])
    return pop_snvs


def parse_snp_id(snp_id):
    items = snp_id.split('||')
    contig = items[0] 
    position = items[1]
    ref = items[2]
    alt = items[3]
    return contig, position, ref, alt

def convert_snv_catalog_to_vcf(species, data_batch):
    vcf_base = os.path.join(config.vcf_path, data_batch)
    if not os.path.exists(vcf_base):
        os.makedirs(vcf_base)
    input_file = os.path.join(config.snv_catalog_path, data_batch, species, 'output', '{}.catalog.noAuto.wtRef.tsv'.format(species))
    output_file = os.path.join(vcf_base, '{}.vcf'.format(species))

    with open(input_file, 'r') as f, open(output_file, 'w') as fout:
        print(f"Converting {input_file} to VCF format")
        header = f.readline().strip().split('\t')
        # Write the VCF header
        fout.write("##fileformat=VCFv4.2\n")
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(header[1:]) + "\n")
        
        line_count = 0
        for line in f:
            parts = line.strip().split('\t')
            snp_info = parse_snp_id(parts[0])
            if snp_info:
                contig, position, ref, alt = snp_info
                # Prepare the genotype info
                genotypes = []
                for gt in parts[1:]:
                    if gt == "255":
                        genotypes.append(".")
                    else:
                        genotype = f"{gt}"
                        genotypes.append(genotype)
                # Write the VCF entry
                fout.write(f"{contig}\t{position}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(genotypes) + "\n")
                line_count += 1
    print(f"VCF file written to {output_file}")
    print("Write {} snps".format(line_count))