import numpy as np
import random
import logging

import config
from utils import snv_utils, metadata_utils

def window_prep(pairs_to_process, snp_file_path, cover_file_path, get_snv_vectors, baselength=40):
    """
    Prepare the input files for combo_prep.py from the output of get_snp_vector
    Input:
    pairs_to_process: list of pairs to process
    snp_file_path: path to the output file with SNPs
    cover_file_path: path to the output file with coverage
    pair_to_snvs: Function that returns a tuple of two vectors: the coverage vector, and the location of SNVs along the coverage vector
    baselength: the length of the window for coverage

    """
    # output: the combo_prep.py format in MAGIC
    snp_file = open(snp_file_path, 'w')
    coverage = open(cover_file_path, 'w')
    
    total_snvs = 0
    total_genome_len = 0
    total_covered_len = 0
    
    cumu_pos = 0 # the position along the concatenated genome
    nr_called = 0
    nr_called_window = 0

    for pair in pairs_to_process:
        pair_snv_locs, pair_coverage = get_snv_vectors(pair[0], pair[1])
        assert(pair_coverage[pair_snv_locs].all())
        curr_pos = 0

        logging.info("Total {} snvs in current pair; genome length {}".format(len(pair_snv_locs), len(pair_coverage)))
        total_snvs += len(pair_snv_locs)
        total_genome_len += len(pair_coverage)
        total_covered_len += pair_coverage.sum()
        
        while curr_pos < len(pair_coverage):
            curr_pos += 1
            cumu_pos += 1
            
            # in combo_prep, position is 1 based
            if pair_coverage[curr_pos-1]:
                # nr_called is the number covered sites in this pair between snvs
                # nr_called_window is the number of covered sites in the current window of length baselength
                nr_called += 1
                nr_called_window += 1
            if cumu_pos % baselength == 0:
                # save current window to coverfile
                # since using cumu_pos, windowing is done with the concatenated genomes, not separately for each pair
                
                # since python2, need other syntax...
                coverage.write(str(nr_called_window)+"\n")
#                 print(nr_called_window, file=coverage)
                nr_called_window = 0
            if (curr_pos-1) in pair_snv_locs:
                # the current position is a segregating site
                # here for QP data, we just put a filler for the allele field
                alleles = '01'
                chromo = 'NA'
                # the snp position must be cumu_pos for windower to work
                snp_file.write("\t".join(map(str, [chromo, cumu_pos, nr_called, alleles])) + "\n")
#                 print(chromo, cumu_pos, nr_called, alleles, sep="\t", file=snp_file)
                # start counting the sites in the next interval between SNVs
                nr_called = 0
        
        # note: there could be windows after the last SNV; they will be empty
        # combo_prep.py only reads to the end of the window that contains the last SNV
    logging.info("Found {} SNVs; total genome length: {}; total covered length: {}".format(total_snvs, total_genome_len, total_covered_len))
    snp_file.close()
    coverage.close()


def sample_random_pairs(indices, num_pairs):
    """
    Randomly sample pairs of indices from the given set, ensuring no index is reused.
    Each index is used only once and each pair is unique.

    Parameters:
        indices (set or list): A set or list of indices from which to sample pairs.
        num_pairs (int): The number of pairs to sample.

    Returns:
        list of tuples: A list containing the randomly sampled pairs.
        None if the number of requested pairs is too large.

    Raises:
        ValueError: If the number of requested pairs is not feasible with the given indices.
    """
    # Convert indices to a list and shuffle it for randomness
    shuffled_indices = list(indices)
    random.shuffle(shuffled_indices)

    # Check if the number of requested pairs is feasible
    if num_pairs > len(shuffled_indices) // 2:
        raise ValueError("The number of requested pairs is too large for the provided set of indices.")

    # Generate the pairs from the shuffled list
    pairs = [(shuffled_indices[i], shuffled_indices[i+1]) for i in range(0, 2*num_pairs, 2)]
    return pairs


def sample_cross_set_pairs(set1, set2, num_pairs):
    """
    Randomly sample pairs of indices such that each pair contains one index from each of the two sets.
    Each index is used only once and each pair is unique.

    Parameters:
        set1 (set or list): A set or list of indices for the first element of each pair.
        set2 (set or list): A set or list of indices for the second element of each pair.
        num_pairs (int): The number of pairs to sample.

    Returns:
        list of tuples: A list containing the randomly sampled pairs.
        None if the number of requested pairs is too large.

    Raises:
        ValueError: If the number of requested pairs is not feasible with the given indices.
    """
    # Shuffle both lists independently for randomness
    list1 = list(set1)
    list2 = list(set2)
    random.shuffle(list1)
    random.shuffle(list2)

    # Check if the number of requested pairs is feasible
    if num_pairs > min(len(list1), len(list2)):
        raise ValueError("The number of requested pairs exceeds the size of the smaller set.")

    # Generate the pairs from the shuffled lists
    pairs = [(list1[i], list2[i]) for i in range(num_pairs)]
    return pairs