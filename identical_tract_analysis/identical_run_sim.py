"""
Functions for simulating how identical run distribution changes 
with mutation accumulation
"""

import numpy as np
import pandas as pd

def simulate_segments(L, mutation_rate):
    """
    For each segment of length L, simulate how mutations break it into 
    smaller segments

    Parameters
    ----------
    L : int
        Length of the segment
    mutation_rate : float
        Scaled mutation rate per site
    """
    # Step 1: Generate the number of mutations N using Poisson distribution
    N = np.random.poisson(mutation_rate * L)
    if N == 0:
        # If there are no mutations, the longest segment is the entire genome
        return [L]
    # Step 2: Randomly place N mutations along the segment of length L
    mutation_positions = np.sort(np.random.randint(0, L, N))
    # Step 3: Compute the gaps between consecutive mutations
    # Include the segments before the first and after the last mutation
    gaps = np.diff(np.concatenate(([0], mutation_positions, [L])))
    return gaps

def simulate_mut_accumulation_one_genome(Ls, mutation_rate):
    """
    Take the identical segments for a simple pairwise comparison
    and simulate how they change with mutation accumulation
    """
    segments_lengths = [simulate_segments(L, mutation_rate) for L in Ls]
    return np.concatenate(segments_lengths)