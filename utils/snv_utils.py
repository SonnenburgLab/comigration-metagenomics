"""
Utilities for manipulating SNV catalogs.
"""

import os
import pandas as pd
import numpy as np
import logging
from pathlib import Path

import config
from utils import metadata_utils

class SNVHelper:
    def __init__(self, species_name, data_batch, metadata_helper=None, cache_snvs=True):
        self.species_name = species_name
        self.data_batch = data_batch
        self.snv_file = config.snv_catalog_path / self.data_batch / self.species_name / 'output' / f'{self.species_name}.catalog.noAuto.wtRef.tsv'
        self.degeneracy_file = config.snv_catalog_path / self.data_batch / self.species_name / 'output' / f'{self.species_name}_4D_sites.tsv'
        self.core_gene_file = config.snv_catalog_path / self.data_batch / self.species_name / 'output' / f'{self.species_name}_core_genome_mask.tsv'
        if metadata_helper is None:
            self.metahelper = metadata_utils.MetadataHelper(data_batch)
        else:
            assert data_batch == metadata_helper.data_batch

        logging.info(f"Loading SNV data for {species_name}")
        if cache_snvs:
            self.cache_format = config.cache_format
            self.cache_syn_snvs_path = config.cache_snvs_path / data_batch / 'syn_snvs' / f'{species_name}_syn_snvs.{self.cache_format}'
            # check if already computed
            if self.cache_syn_snvs_path.exists():
                logging.info(f"Loading cached SNVs from {self.cache_syn_snvs_path}")
                self.syn_snvs = self.load_cached_snvs(self.cache_syn_snvs_path, self.cache_format)
                self.check_reference_syn_snvs()
                self.genome_names = self.syn_snvs.columns.values.tolist()
                return
    
        logging.info("Computing synonymous SNVs")
        self.load_degeneracy()
        self.load_core_gene_index()
        # this is usually the slower step because it involves parsing a big csv
        self.load_snvs()
        self.syn_snvs = self.compute_4D_core_snvs()
        self.check_reference_syn_snvs()
        if cache_snvs:
            self.cache_snvs()
        self.genome_names = self.syn_snvs.columns.values.tolist()

    def load_degeneracy(self):
        # degeneracy df only contains sites at the last codon position, so does not contain all sites
        self.degeneracy = load_degeneracy(self.degeneracy_file)
    
    def load_snvs(self):
        self.snvs = load_snvs(self.snv_file)

    def load_core_gene_index(self):
        self.core_gene_index = load_core_gene_mask(self.core_gene_file)

    def load_cached_coverage_matrix(self):
        self.coverage_matrix_path = config.cache_snvs_path / self.data_batch / 'coverage' / f'{self.species_name}_coverage.{self.cache_format}'
        if not self.coverage_matrix_path.exists():
            raise RuntimeError(f'Coverage matrix not found at {self.coverage_matrix_path}')
        elif self.cache_format == 'feather':
            self.coverage = pd.read_feather(self.coverage_matrix_path)
        elif self.cache_format == 'parquet':
            self.coverage = pd.read_parquet(self.coverage_matrix_path)
        self.coverage.set_index(['Contig', 'Site'], inplace=True)
        self.coverage.sort_index(inplace=True)
    
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
        cache_syn_snvs_path = config.cache_snvs_path / self.data_batch / 'syn_snvs'
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

    def get_core_gene_index(self):
        # NOT sorted!
        if not hasattr(self, 'core_gene_index'):
            # loading necessary data
            logging.info('Data not loaded; loading raw data')
            self.load_core_gene_index()
        return self.core_gene_index
    
    def get_degeneracy(self):
        if not hasattr(self, 'degeneracy'):
            # loading necessary data
            logging.info('Data not loaded; loading raw data')
            self.load_degeneracy()
        return self.degeneracy
    
    def get_coverage(self):
        if not hasattr(self, 'coverage'):
            # loading necessary data
            logging.info('Data not loaded; loading cached coverage matrix')
            self.load_cached_coverage_matrix()
        return self.coverage

    def get_4D_core_indices(self):
        # NOT sorted!
        #TODO: using multiindexing; could be a bit slow compared to boolean masks
        core_index = self.get_core_gene_index()
        if core_index is None:
            return None
        degen = self.get_degeneracy()
        degen_4D = degen[degen['Degeneracy'] == 4]
        core_mask = degen_4D.index.isin(core_index)
        return degen_4D[core_mask].index
    
    def get_4D_core_genome_length(self):
        core_index = self.get_4D_core_indices()
        if core_index is None:
            return 0
        return len(core_index)
    
    def get_num_syn_snps(self):
        return self.syn_snvs.shape[0]
    
    def compute_4D_core_snvs(self):
        """
        The resulting DataFrame has the following multi-index, which specify a unique snv:
        contig, position, ref, alt

        The values of the DataFrame should be 0, 1, or 255, where 0 is the reference allele, 1 is the alternate allele, and 255 is missing data.
        """
        if not hasattr(self, 'snvs'):
            # loading necessary data
            logging.info('Data not loaded; loading raw data')
            self.load_snvs()

        core_indices = self.get_4D_core_indices()
        if core_indices is None:
            raise ValueError('No core genome found')
        snv_mask = self.snvs.index.isin(core_indices)
        syn_snvs = self.snvs[snv_mask].copy()
        syn_snvs.set_index(['Ref', 'Alt'], append=True, inplace=True)
        syn_snvs.drop(['snp_id'], axis=1, inplace=True)
        return syn_snvs.sort_index()
    
    def check_reference_syn_snvs(self):
        # check if the reference genome is in syn_snvs, and if so, make sure it's set to all 0s
        # the reason this function exists is that databatch 240509 includes the reference in the snv catalog as a column of 255s
        rep_genome = self.metahelper.get_species_rep_genome(self.species_name)
        if rep_genome in self.syn_snvs.columns:
            self.syn_snvs[rep_genome] = 0
    
    def get_biallelic_core_snvs(self):
        # TODO: implement by grouping the snvs by index and checking if there are only two alleles
        pass

    def get_population_mask(self, pop_name, allowed_mags=None):
        # need MetadataHelper
        self.genome_names = self.syn_snvs.columns.values.tolist()
        if allowed_mags is not None:
            pop_mask = np.array([(self.metahelper.get_mag_pop(x) == pop_name) and (x in allowed_mags) for x in self.genome_names])
        else:
            pop_mask = np.array([self.metahelper.get_mag_pop(x) == pop_name for x in self.genome_names])
        return pop_mask
    
    def get_population_snvs(self, pop_name):
        pop_mask = self.get_population_mask(pop_name)
        return self.syn_snvs.loc[:, pop_mask]
    
    def compute_population_coverage(self, pop_name, allowed_mags=None):
        # default: computing for 4D core sites
        # return: DataFrame of num reference & alternate alleles for each snv
        pop_mask = self.get_population_mask(pop_name, allowed_mags)
        if pop_mask.sum() == 0:
            res = pd.DataFrame(index=self.syn_snvs.index)
            res['Ref count'] = 0
            res['Alt count'] = 0
            return res
        pop_snvs = self.syn_snvs.loc[:, pop_mask]
        ref_counts = (pop_snvs==0).sum(axis=1)
        alt_counts = (pop_snvs==1).sum(axis=1)
        return pd.DataFrame({'Ref count': ref_counts, 'Alt count': alt_counts}, index=pop_snvs.index)

    def save_dadi_data_dict(self, pops, allowed_mags=None, output_path=None):
        """
        Write the snps to a file in the format required by dadi
        # https://dadi.readthedocs.io/en/latest/user-guide/importing-data/

        """
        if output_path is None:
            # output file default path
            output = config.sfs_path / self.data_batch
        else:
            output = Path(output_path)
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
        # these last ones are joined together as the snp_id
        header_items = header_items + ['Contig', 'Position', 'snp_id_Ref', 'snp_id_Alt', '\n']

        # calculating counts for snvs in each population
        snv_counts = {}
        for pop in pops:
            snv_counts[pop] = self.compute_population_coverage(pop, allowed_mags)
        
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
                
                line_items = [ref_string, out_string, a1] + refs + [a2] + alts + [contig, pos, a1, a2, '\n']
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


# a child class for handling SNV vectors when doing between MAG comparison
class PairwiseSNVHelper(SNVHelper):
    """
    A child class of SNVHelper that is specifically designed for handling SNV vectors.

    Loads additional data, such as the coverage matrix, and provides additional methods for computing SNV vectors.
    """
    def __init__(self, species_name, data_batch):
        super().__init__(species_name, data_batch, metadata_helper=None, cache_snvs=True)
        self.load_cached_coverage_matrix()
        core_4D_sites = self.get_4D_core_indices()
        core_to_4D = self.coverage.index.isin(core_4D_sites)
        self.core_4D_coverage = self.coverage[core_to_4D].copy()

        self.compute_bi_allelic_snvs()

    def compute_bi_allelic_snvs(self):
        """
        Find the bi-allelic SNVs in the SNV catalog and save the locations of them.
        Mask the multi-allelic sites.
        """
        # iterate over the first two levels of syn_snvs index
        multi_sites = self.syn_snvs.groupby(level=[0, 1]).filter(lambda x: len(x) > 1).index
        multi_sites_index = multi_sites.droplevel(level=[2, 3]).drop_duplicates()
        snv_to_bi_mask = ~(self.syn_snvs.index.isin(multi_sites))
        self.bi_snvs = self.syn_snvs[snv_to_bi_mask]

        multi_mask = self.core_4D_coverage.index.isin(multi_sites_index)
        print("{} multi-allelic sites masked".format(multi_mask.sum()))
        self.core_4D_coverage.loc[multi_mask] = 0 # masking the multi-allelic sites for now; TODO revisit later

        self.bi_snv_index = self.bi_snvs.index.droplevel(level=[2, 3])
        self.core_4D_to_bi_snvs = self.core_4D_coverage.index.isin(self.bi_snv_index)

        self.bi_snv_vector_locs = np.nonzero(self.core_4D_to_bi_snvs)[0]

    def get_snv_vector(self, sample1, sample2):
        """
        Compute the location of SNVs that are different between two samples.
        """
        cover_mat = self.core_4D_coverage[[sample1, sample2]].values
        cover_vec = np.all(cover_mat==1, axis=1)

        snv_cover = cover_vec[self.core_4D_to_bi_snvs]

        is_snv = ((self.bi_snvs[sample1]==0) & (self.bi_snvs[sample2]==1)) | ((self.bi_snvs[sample1]==1) & (self.bi_snvs[sample2]==0))
        is_snv = is_snv.values & snv_cover 
        snv_locs = self.bi_snv_vector_locs[is_snv].copy()
        return snv_locs, cover_vec

    @staticmethod
    def compute_identical_fraction(snp_vec, block_size=1000):
        snp_blocks = to_snv_blocks(snp_vec, block_size)
        frac_id = np.mean(snp_blocks == 0)
        return frac_id
    
    @staticmethod
    def compute_runs(snp_vec, positions=None, return_locs=False):
        """
        Compute the runs of identical sites given a vector of 0s and 1s.
        Runs include start->first snp and last snp->end
        If coordinate position of each site is provided, the function will compute runs in terms of coordinate positions.
        If return_locs is True, the function will return the start and end locations of each run.
        """
        # get the locations of snps in the vector
        padded_vec = np.ones(len(snp_vec) + 2)
        padded_vec[1:-1] = snp_vec
        # the locations here are along the snp_vec, not the reference genome
        site_locations = np.nonzero(padded_vec)[0]
        if positions is not None:
            padded_locs = np.zeros(len(positions) + 2)
            padded_locs[0] = positions[0] - 1
            padded_locs[1:-1] = positions 
            padded_locs[-1] = positions[-1] + 1
            locs = padded_locs[site_locations]
        else:
            locs = site_locations
        runs = site_locations[1:] - site_locations[:-1] - 1
        starts = locs[:-1][runs > 0]
        ends = locs[1:][runs > 0]
        runs = runs[runs > 0]
        if return_locs:
            return runs, starts, ends
        return runs
    
    @staticmethod
    def compute_runs_by_contigs(snp_vec, contigs, positions=None, return_locs=False):
        """
        Compute the runs of identical sites given a vector of 0s and 1s, separated by contigs.
        That is, runs are computed separately for each contig, and no run can span two contigs.
        Assuming that snp_vec is ordered by contigs.

        Parameters:
        - snp_vec: numpy array of 0s and 1s
        - contigs: numpy array of contig names, same shape as snp_vec
        - positions: numpy array of positions along the reference genome, same shape as snp_vec
        - return_locs: if True, return the start and end locations of each run
        """
        all_runs = []
        all_starts = []
        all_ends = []
        for contig in pd.unique(contigs):
            subvec = snp_vec[contigs==contig]
            if positions is not None:
                subloc = positions[contigs==contig]
            else:
                subloc = None
            res = PairwiseSNVHelper.compute_runs(subvec, subloc, return_locs=True)
            all_runs.append(res[0])
            all_starts.append(res[1])
            all_ends.append(res[2])
        if return_locs:
            return np.concatenate(all_runs), np.concatenate(all_starts), np.concatenate(all_ends)
        else:
            return np.concatenate(all_runs)
        
    def compute_runs_for_pairs(self, pairs, separate_contigs=False):
        """
        Compute the runs of identical sites for each pair of samples / MAGs / genomes.
        If separate_contigs is True, the runs are computed separately for each contig.
        """
        if separate_contigs:
            raise NotImplementedError('Separate contigs not implemented yet')
        runs = {}
        for pair in pairs:
            snv_locs, cover_vec = self.get_snv_vector(pair[0], pair[1])
            snv_vec = np.zeros(cover_vec.shape)
            snv_vec[snv_locs] = 1
            snv_vec = snv_vec[cover_vec]
            runs[pair] = self.compute_runs(snv_vec)
        return runs
    
    def compute_frac_id_for_pairs(self, pairs, block_size=1000):
        """
        Compute the fraction of identical sites for each pair of samples / MAGs / genomes.
        """
        genome1s = [pair[0] for pair in pairs]
        genome2s = [pair[1] for pair in pairs]
        res = pd.DataFrame()
        res['genome1'] = genome1s
        res['genome2'] = genome2s
        res['frac_id'] = np.nan
        for i, pair in enumerate(pairs):
            snv_locs, cover_vec = self.get_snv_vector(pair[0], pair[1])
            snv_vec = np.zeros(cover_vec.shape)
            snv_vec[snv_locs] = 1
            snv_vec = snv_vec[cover_vec]
            frac_id = self.compute_identical_fraction(snv_vec, block_size)
            res.at[i, 'frac_id'] = frac_id
        return res
        
    @staticmethod
    def generate_pairs(samples):
        # helper for doing pairwise comparisons
        # generate all unique pairs within a list
        pairs = []
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                pairs.append((samples[i], samples[j]))
        return pairs

    @staticmethod
    def generate_pairs_between(samples1, samples2):
        # generate all unique pairs between two lists
        pairs = []
        for i in range(len(samples1)):
            for j in range(len(samples2)):
                pairs.append((samples1[i], samples2[j]))
        return pairs
    
    def get_all_genome_pairs(self):
        return self.generate_pairs(self.genome_names)


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

def load_core_gene_mask(core_gene_file):
    # check if file is empty
    try:
        core_gene_df = pd.read_csv(core_gene_file, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        return None
    core_gene_df.columns = ['Contig', 'Site']
    core_gene_index = core_gene_df.set_index(['Contig', 'Site']).index.drop_duplicates()
    return core_gene_index

# next few functions are for handling coords files

def filename_to_samples(filename):
    items = filename.stem.split('-')
    return items[0], items[1]

def parse_coords(filename):
    """
    Parse the coordinates from the output of nucmer (which is part of the UHGG snv calling pipeline).
    Example first few lines:

    Agathobacter_rectalis__TS_ADULT_82.fa Agathobacter_rectalis/genomes/Agathobacter_rectalis__ERS396293.fa
    NUCMER

        [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
    =====================================================================================
    59522    61379  |        1     1858  |     1858     1858  |    98.98  | TS_ADULT_82__NODE_593_length_71722_cov_112.888729	ERS396293|6|k85_102948

    """
    data = []
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if i < 5:  # Skip the header lines
                continue
            # Split the line on '|' first, then on whitespace and tab as appropriate
            parts = line.split('|')
            ref_coord = parts[0].split()
            query_coord = parts[1].split()
            lengths = parts[2].split()
            percent_identity = parts[3].strip()

            # now rejoin the tags, since they may contain '|' characters
            tags = '|'.join(parts[4:])
            items = tags.split()
            ref_contig = items[0]
            query_contig = items[1]
            row = ref_coord + query_coord + lengths + [percent_identity, ref_contig, query_contig]
            data.append(row)
    df = pd.DataFrame(data, columns=['ref_start', 'ref_end', 'query_start', 'query_end', 'ref_length', 'query_length', 'percent_identity', 'ref_contig', 'query_contig'])
    df[['ref_start', 'ref_end', 'query_start', 'query_end', 'ref_length', 'query_length']] = df[['ref_start', 'ref_end', 'query_start', 'query_end', 'ref_length', 'query_length']].astype(int)
    df['percent_identity'] = df['percent_identity'].astype(float)
    return df

def compute_coverage_matrix_from_coords(filenames, core_gene_idx, mag_names):
    """
    Compute the coverage matrix from the coordinates files output by nucmer.
    Input:
    - filenames: list of filenames for all the coordinates files; name should be in the format ref_name-query_name.coords
    - core_gene_idx: index of the core genes, as a pandas MultiIndex
    - mag_names: list of MAG names
    Output:
    - coverage_matrix: DataFrame of shape (len(core_gene_idx), len(mag_names)) with the coverage for each MAG at each core gene
    Values are integers, representing the number of MAG chunks covering a given reference location.
    """
    # just in case; need to make sure sites are ordered and unique
    core_gene_idx = core_gene_idx.drop_duplicates()
    core_gene_idx = core_gene_idx.sort_values()
    mask = np.empty((len(core_gene_idx), len(mag_names)), dtype=int)
    coverage_matrix = pd.DataFrame(mask, index=core_gene_idx, columns=mag_names)

    # saving time by caching the contig masks
    contigs = coverage_matrix.index.get_level_values(0)
    core_contigs = contigs.unique()
    contig_to_mask = {contig: contigs == contig for contig in core_contigs}

    sites = coverage_matrix.index.get_level_values(1)

    for filename in filenames:
        logging.info(f'Processing {filename}')

        ref_name, query_name = filename_to_samples(filename)
        if query_name not in mag_names:
            logging.warning(f'Genome {query_name} not found in provided mag_names')
            continue
        coords_df = parse_coords(filename)

        # saves panda indexing overhead by setting a numpy array first, and finally assigning it to the coverage matrix
        coverage = np.zeros(len(core_gene_idx), dtype=int)
        for contig, chunks in coords_df.groupby('ref_contig'):
            if contig not in core_contigs:
                continue
            contig_mask = contig_to_mask[contig]
            for i, row in chunks.iterrows():
                ref_start = row['ref_start']
                ref_end = row['ref_end']
                chunk_mask = contig_mask & (sites >= ref_start) & (sites <= ref_end)
                coverage[chunk_mask] += 1
        coverage_matrix.loc[:, query_name] = coverage
    return coverage_matrix


# A few functions copied from Liu & Good 2024 microbiome recombination project
def length_to_num_blocks(seq_len, block_size):
    # Magical formula that works for all cases
    return (seq_len + block_size - 1) // block_size

def to_snv_blocks(bool_array, block_size):
    """
    Converting a boolean array into blocks of True counts. The last
    block could be shorter than block_size
    :param bool_array:
    :param block_size:
    :return: An array of counts of Trues in blocks
    """
    # coarse-graining the bool array (snp array) into blocks
    num_blocks = length_to_num_blocks(len(bool_array), block_size)
    bins = np.arange(0, num_blocks * block_size + 1, block_size)
    counts, _ = np.histogram(np.nonzero(bool_array), bins)
    return counts