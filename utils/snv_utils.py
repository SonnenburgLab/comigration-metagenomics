"""
Utilities for manipulating SNV catalogs.
"""

import os
import pandas as pd
import numpy as np
import logging

import config
from utils import metadata_utils

class SNVHelper:
    def __init__(self, species_name, data_batch, metadata_helper=None, cache_snvs=False):
        self.species_name = species_name
        self.data_batch = data_batch
        if metadata_helper is None:
            self.metahelper = metadata_utils.MetadataHelper(data_batch)
        else:
            assert data_batch == metadata_helper.data_batch

        logging.info(f"Loading SNV data for {species_name}")
        if cache_snvs:
            self.cache_format = config.cache_format
            cache_syn_snvs_path = config.cache_snvs_path / data_batch / f'{species_name}_syn_snvs.{self.cache_format}'
            # check if already computed
            if cache_syn_snvs_path.exists():
                logging.info(f"Loading cached SNVs from {cache_syn_snvs_path}")
                self.syn_snvs = self.load_cached_snvs(cache_syn_snvs_path, self.cache_format)
                self.genome_names = self.syn_snvs.columns.values.tolist()
                return
        logging.info("Computing synonymous SNVs")
        self.load_raw_snvs()
        self.syn_snvs = self.get_4D_core_snvs()
        if cache_snvs:
            self.cache_snvs()
        self.genome_names = self.syn_snvs.columns.values.tolist()

    def load_raw_snvs(self):
        snv_file = config.snv_catalog_path / self.data_batch / self.species_name / 'output' / f'{self.species_name}.catalog.noAuto.wtRef.tsv'
        degeneracy_file = config.snv_catalog_path / self.data_batch / self.species_name / 'output' / f'{self.species_name}_4D_sites.tsv'
        core_gene_file = config.snv_catalog_path / self.data_batch / self.species_name / 'output' / f'{self.species_name}_core_genome_mask.tsv'

        # TODO: move these to static methods?
        self.snvs = load_snvs(snv_file)
        # degeneracy df only contains sites at the last codon position, so does not contain all sites
        self.degeneracy = load_degeneracy(degeneracy_file)
        self.core_gene_index = load_core_gene_mask(core_gene_file)
    
    # static helper for loading cached dataframes
    @staticmethod
    def load_cached_snvs(path, format):
        if format not in ['feather', 'parquet']:
            raise ValueError('format should be either feather or parquet')
        if format=='feather':
            df = pd.read_feather(path)
        elif format=='parquet':
            df = pd.read_feather(path)
        # the following code is not necessary if the data is saved correctly
        # for col in df.columns:
        #     if df[col].dtype == 'O':
        #         if isinstance(df[col][0], str):
        #             continue
        #         elif isinstance(df[col][0], bytes):
        #             # some of string series (from python2) were saved as bytes
        #             df[col] = df[col].str.decode('utf-8')
        df.set_index(['Contig', 'Site', 'Ref', 'Alt'], inplace=True)
        return df
    
    def cache_snvs(self):
        cache_syn_snvs_path = config.cache_snvs_path / self.data_batch 
        if not cache_syn_snvs_path.exists():
            os.makedirs(cache_syn_snvs_path)
        cache_syn_snvs_path = cache_syn_snvs_path / f'{self.species_name}_syn_snvs.{self.cache_format}'
        logging.info(f"Caching SNVs to {cache_syn_snvs_path}")
        if self.cache_format not in ['feather', 'parquet']:
            raise ValueError('format should be either feather or parquet')
        if self.cache_format=='feather':
            self.syn_snvs.reset_index().to_feather(cache_syn_snvs_path)
        elif self.cache_format=='parquet':
            self.syn_snvs.reset_index().to_parquet(cache_syn_snvs_path)

    def get_4D_core_indices(self):
        #TODO: using multiindexing; could be a bit slow compared to boolean masks
        if self.core_gene_index is None:
            return None
        degen_4D = self.degeneracy[self.degeneracy['Degeneracy'] == 4]
        core_mask = degen_4D.index.isin(self.core_gene_index)
        return degen_4D[core_mask].index
    
    def get_4D_core_snvs(self):
        """
        The resulting DataFrame has the following multi-index, which specify a unique snv:
        contig, position, ref, alt

        The values of the DataFrame should be 0, 1, or 255, where 0 is the reference allele, 1 is the alternate allele, and 255 is missing data.
        """
        # check if core_gene_index is in attributes

        if not hasattr(self, 'core_gene_index'):
            raise ValueError('Raw data not loaded; call load_raw_snvs() first')
        core_indices = self.get_4D_core_indices()
        if core_indices is None:
            raise ValueError('No core genome found')
        snv_mask = self.snvs.index.isin(core_indices)
        syn_snvs = self.snvs[snv_mask].copy()
        syn_snvs.set_index(['Ref', 'Alt'], append=True, inplace=True)
        syn_snvs.drop(['snp_id'], axis=1, inplace=True)
        return syn_snvs.sort_index()
    
    def get_biallelic_core_snvs(self):
        # TODO: implement by grouping the snvs by index and checking if there are only two alleles
        pass

    def get_population_mask(self, pop_name):
        # need MetadataHelper
        pops = np.array([self.metahelper.get_mag_pop(x) for x in self.genome_names])
        return pops == pop_name
    
    def compute_population_coverage(self, pop_name):
        # default: computing for 4D core sites
        # return: DataFrame of num reference & alternate alleles for each snv
        pop_mask = self.get_population_mask(pop_name)
        if pop_mask.sum() == 0:
            res = pd.DataFrame(index=self.syn_snvs.index)
            res['Ref count'] = 0
            res['Alt count'] = 0
            return res
        pop_snvs = self.syn_snvs.loc[:, pop_mask]
        ref_counts = (pop_snvs==0).sum(axis=1)
        alt_counts = (pop_snvs==1).sum(axis=1)
        return pd.DataFrame({'Ref count': ref_counts, 'Alt count': alt_counts}, index=pop_snvs.index)

    def save_dadi_data_dict(self, pops):
        """
        Write the snps to a file in the format required by dadi
        # https://dadi.readthedocs.io/en/latest/user-guide/importing-data/

        """
        # output file default path
        output = config.sfs_path / self.data_batch
        if not os.path.exists(output):
            os.makedirs(output)
        output_file = output / f'{self.species_name}.snps.txt'

        # can arbitrary number of populations; preparing header
        header_items = ['Ref', 'Out', 'Allele1']
        for pop in pops:
            header_items.append(pop)
        header_items.append('Allele2')
        for pop in pops:
            header_items.append(pop)
        header_items = header_items + ['Contig', 'Position', '\n']

        # calculating counts for snvs in each population
        snv_counts = {}
        for pop in pops:
            snv_counts[pop] = self.compute_population_coverage(pop)
        
        # Save to file; iterating over snvs
        logging.info(f'Writing to file: {output_file}')
        with open(output_file, 'w') as snp_file:
            snp_file.write('\t'.join(header_items))

            for ind in self.syn_snvs.index:
                # iterate over snp rows
                # first collect snp information
                ref_string = '-'+ind[2]+'-'
                out_string = '---'
                # ref and alt alleles
                a1 = ind[2]
                a2 = ind[3]
                # contig and position; only useful for tagging the snv
                pos = ind[1]
                contig = ind[0]

                refs = []
                alts = []
                # now write over populations
                for pop in pops:
                    refs.append(snv_counts[pop].at[ind, 'Ref count'])
                    alts.append(snv_counts[pop].at[ind, 'Alt count'])
                
                line_items = [ref_string, out_string, a1] + refs + [a2] + alts + [contig, pos, '\n']
                line_items = [str(x) for x in line_items]
                snp_file.write('\t'.join(line_items))
        logging.info(f'File written to {output_file}')

    def save_vcf(self):
        """
        Save the 4D core sites to a VCF file.
        """
        vcf_base = config.vcf_path / self.data_batch
        if not vcf_base.exists():
            vcf_base.mkdir()
        output_file = vcf_base / f'{self.species_name}.vcf'

        with open(output_file, 'w') as fout:
            logging.info(f"Converting {self.syn_snvs.shape[0]} snvs to VCF format")
            # Write the VCF header
            fout.write("##fileformat=VCFv4.2\n")
            fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(self.syn_snvs.columns) + "\n")
            
            for ind, line in self.syn_snvs.iterrows():
                contig, position, ref, alt = ind
                genotypes = []
                for gt in line:
                    if gt == 255:
                        genotypes.append(".")
                    else:
                        genotype = f"{gt}"
                        genotypes.append(genotype)
                fout.write(f"{contig}\t{position}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(genotypes) + "\n")
        logging.info(f"VCF file written to {output_file}")

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
    snv_catalog['Ref'] = snv_catalog['snp_id'].apply(lambda x: x.split('||')[-2])
    snv_catalog['Alt'] = snv_catalog['snp_id'].apply(lambda x: x.split('||')[-1])
    snv_catalog.set_index(['Contig', 'Site'], inplace=True)
    return snv_catalog

def convert_snv_catalog_to_vcf(species, data_batch, syn_core_only=True):
    vcf_base = os.path.join(config.vcf_path, data_batch)
    if not os.path.exists(vcf_base):
        os.makedirs(vcf_base)
    input_file = os.path.join(config.snv_catalog_path, data_batch, species, 'output', '{}.catalog.noAuto.wtRef.tsv'.format(species))
    output_file = os.path.join(vcf_base, '{}.vcf'.format(species))
    if syn_core_only:
        core_gene_file = os.path.join(config.snv_catalog_path, data_batch, species, 'output', '{}_core_genome_mask.tsv'.format(species))
        degeneracy_file = os.path.join(config.snv_catalog_path, data_batch, species, 'output', '{}_4D_sites.tsv'.format(species))
        snv_catalog, _  = load_and_filter_snv_catalog(input_file, degeneracy_file, core_gene_file)
    else:
        snv_catalog = load_snvs(input_file)

    with open(input_file, 'r') as f, open(output_file, 'w') as fout:
        print(f"Converting {input_file} to VCF format")
        header = f.readline().strip().split('\t')
        # Write the VCF header
        fout.write("##fileformat=VCFv4.2\n")
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(header[1:]) + "\n")
        
        line_count = 0
        for ind, line in snv_catalog.iterrows():
            contig, position, ref, alt = ind
            genotypes = []
            for gt in line[header[1:]]:
                if gt == 255:
                    genotypes.append(".")
                else:
                    genotype = f"{gt}"
                    genotypes.append(genotype)
            fout.write(f"{contig}\t{position}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(genotypes) + "\n")
            line_count += 1
    print(f"VCF file written to {output_file}")
    print("Write {} snps".format(line_count))

def load_core_gene_mask(core_gene_file):
    # check if file is empty
    try:
        core_gene_df = pd.read_csv(core_gene_file, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        return None
    core_gene_df.columns = ['Contig', 'Site']
    core_gene_index = core_gene_df.set_index(['Contig', 'Site']).index
    return core_gene_index

def load_dadi_snps_metadata(batch_name):
    metadata_path = os.path.join(config.dadi_dat_path, '{}__metadata.csv'.format(batch_name))
    return pd.read_csv(metadata_path, index_col=0)

def load_and_filter_snv_catalog(snv_file, degeneracy_file, core_gene_file):
    """
    TODO: deprecated; use SNVHelper instead; also the 1D sites are not full here, only the last codon position
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

    syn_snvs = syn_snvs.sort_index()
    nonsyn_snvs = nonsyn_snvs.sort_index()

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
    TODO: implement a method in SNVHelper
    Return both the total length of the core genome and the length of the 4D core genome.
    """
    core_gene_index = load_core_gene_mask(core_gene_file)
    if core_gene_index is None:
        return 0, 0
    core_4D = get_4D_core_indices(degeneracy_file, core_gene_file)
    return core_gene_index.shape[0], len(core_4D)
