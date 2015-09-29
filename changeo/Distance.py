"""
Distance data structures and functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'
from changeo import __version__, __date__

# Imports
import numpy as np
import pandas as pd
from itertools import combinations, product
from pkg_resources import resource_stream
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

# Presto and changeo imports
from presto.Sequence import scoreDNA, scoreAA

# Load model data
with resource_stream(__name__, 'data/M1N_Distance.tab') as f:
    #m1n_distance = pd.read_csv(f, sep='\t', index_col=0).to_dict()
    m1n_model = pd.read_csv(f, sep='\t', index_col=0)

with resource_stream(__name__, 'data/HS1F_Distance.tab') as f:
    hs1f_model = pd.read_csv(f, sep='\t', index_col=0)

with resource_stream(__name__, 'data/HS5F_Distance.tab') as f:
    hs5f_model = pd.read_csv(f, sep='\t', index_col=0)


def getDNADistMatrix(mat=None, mask_dist=0, gap_dist=0):
    """
    Generates a DNA distance matrix

    Arguments:
      mat : Input distance matrix to extend to full alphabet;
            if unspecified, creates Hamming distance matrix that incorporates
            IUPAC equivalencies
      mask_dist : Distance for all matches against an N character
      gap_dist : Distance for all matches against a gap (-, .) character

    Returns:
      DataFrame : pandas.DataFrame of distances
    """
    IUPAC_chars = list('-.ACGTRYSWKMBDHVN')
    mask_char = 'N'

    # Default matrix to inf
    dist_mat = pd.DataFrame(float('inf'), index=IUPAC_chars, columns=IUPAC_chars,
                            dtype=float)
    # Set gap distance
    for c in '-.':
        dist_mat.loc[c] = dist_mat.loc[:, c] = gap_dist

    # Set mask distance
    dist_mat.loc[mask_char] = dist_mat.loc[:, mask_char] = mask_dist

    # Fill in provided distances from input matrix
    if mat is not None:
        for i,j in product(mat.index, mat.columns):
            dist_mat.at[i, j] = mat.at[i, j]
    # If no input matrix, create IUPAC-defined Hamming distance
    else:
        for i,j in product(dist_mat.index, dist_mat.columns):
            dist_mat.at[i, j] = 1 - scoreDNA(i, j,
                                             mask_score=(1 - mask_dist, 1 - mask_dist),
                                             gap_score=(1 - gap_dist, 1 - gap_dist))

    return dist_mat


def getAADistMatrix(mat=None, mask_dist=0, gap_dist=0):
    """
    Generates an amino acid distance matrix

    Arguments:
      mat : Input distance matrix to extend to full alphabet;
            if unspecified, creates Hamming distance matrix that incorporates
            IUPAC equivalencies
      mask_dict : Score for all matches against an X character
      gap_dist : Score for all matches against a gap (-, .) character

    Returns:
      DataFrame : pandas.DataFrame of distances
    """
    IUPAC_chars = list('-.*ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    mask_char = 'X'

    # Default matrix to inf
    dist_mat = pd.DataFrame(float('inf'), index=IUPAC_chars, columns=IUPAC_chars,
                            dtype=float)

    # Set gap distance
    for c in '-.':
        dist_mat.loc[c] = dist_mat.loc[:, c] = gap_dist

    # Set mask distance
    dist_mat.loc[mask_char] = dist_mat.loc[:, mask_char] = mask_dist

    # Fill in provided distances from input matrix
    if mat is not None:
        for i,j in product(mat.index, mat.columns):
            dist_mat.at[i, j] = mat.at[i, j]
    # If no input matrix, create IUPAC-defined Hamming distance
    else:
        for i,j in product(dist_mat.index, dist_mat.columns):
            dist_mat.at[i, j] = 1 - scoreAA(i, j,
                                            mask_score=(1 - mask_dist, 1 - mask_dist),
                                            gap_score=(1 - gap_dist, 1 - gap_dist))

    return dist_mat


def getNmers(sequences, n):
    """
    Breaks input sequences down into n-mers

    Arguments:
      sequences : List of sequences to be broken into n-mers
      n : Length of n-mers to return

    Returns:
      dict : Dictionary mapping sequence to a list of n-mers
    """
    # Add Ns so first nucleotide is center of first n-mer
    sequences_n = ['N' * ((n - 1) // 2) + seq + 'N' * ((n - 1) // 2) for seq in sequences]
    nmers = {}
    for seq,seqn in zip(sequences,sequences_n):
        nmers[seq] = [seqn[i:i+n] for i in range(len(seqn)-n+1)]
    # nmers = {(seq, [seqn[i:i+n] for i in range(len(seqn)-n+1)]) for seq,seqn in izip(sequences,sequences_n)}

    return nmers


def calcDistances(sequences, n, dist_mat, norm, sym):
    """
    Calculate pairwise distances between input sequences

    Arguments:
      sequences : List of sequences for which to calculate pairwise distances
      n : Length of n-mers to be used in calculating distance
      dist_mat : pandas.DataFrame of mutation distances
      norm : Normalization method
      sym : Symmetry method

    Returns:
      ndarray : numpy matrix of pairwise distances between input sequences
    """
    # Initialize output distance matrix
    dists = np.zeros((len(sequences),len(sequences)))
    # Generate dictionary of n-mers from input sequences
    nmers = getNmers(sequences, n)
    # Iterate over combinations of input sequences
    for j,k in combinations(list(range(len(sequences))), 2):
        # Only consider characters and n-mers with mutations
        mutated = [i for i,(c1,c2) in enumerate(zip(sequences[j],sequences[k])) if c1 != c2]
        seq1 = [sequences[j][i] for i in mutated]
        seq2 = [sequences[k][i] for i in mutated]
        nmer1 = [nmers[sequences[j]][i] for i in mutated]
        nmer2 = [nmers[sequences[k]][i] for i in mutated]

        # Determine normalizing factor
        if norm == 'len':
            norm_by = len(sequences[0])
        elif norm == 'mut':
            norm_by = len(mutated)
        else:
            norm_by = 1

        # Determine symmetry function
        if sym == 'avg':
            sym_fun = np.mean
        elif sym == 'min':
            sym_fun = min
        else:
            sym_fun = sum

        # Calculate distances
        try:
            dists[j, k] = dists[k, j] = \
                    sum([sym_fun([dist_mat.at[c1, n2], dist_mat.at[c2, n1]]) \
                         for c1, c2, n1, n2 in zip(seq1, seq2, nmer1, nmer2)]) / \
                    (norm_by)
        except (KeyError):
            raise KeyError('Unrecognized character in sequence.')

    return dists


def formClusters(dists, link, distance):
    """
    Form clusters based on hierarchical clustering of input distance matrix with
    linkage type and cutoff distance

    Arguments:
      dists : numpy matrix of distances
      link : Linkage type for hierarchical clustering
      distance : Distance at which to cut into clusters

    Returns:
      list : List of cluster assignments
    """
    # Make distance matrix square
    dists = squareform(dists)
    # Compute linkage
    links = linkage(dists, link)
    # Break into clusters based on cutoff
    clusters = fcluster(links, distance, criterion='distance')
    return clusters
