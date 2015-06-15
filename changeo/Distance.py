"""
Distance data structures and functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'
from changeo import __version__, __date__

# Imports
import numpy as np
import pandas as pd
from itertools import combinations, izip, product
from pkg_resources import resource_stream
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

# Presto and changeo imports
from presto.Sequence import scoreDNA, scoreAA

# Data streams
#m1n_data = resource_stream(__name__, 'data/M1N_Distance.tab')
#m3n_data = resource_stream(__name__, 'data/M3N_Distance.tab')
#hs5f_data = resource_stream(__name__, 'data/HS5F_Distance.tab')

# Load model data
with resource_stream(__name__, 'data/M1N_Distance.tab') as f:
    #m1n_distance = pd.read_csv(f, sep='\t', index_col=0).to_dict()
    m1n_distance = pd.read_csv(f, sep='\t', index_col=0)

with resource_stream(__name__, 'data/M3N_Distance.tab') as f:
    m3n_distance = pd.read_csv(f, sep='\t', index_col=0)

with resource_stream(__name__, 'data/HS5F_Distance.tab') as f:
    hs5f_distance = pd.read_csv(f, sep='\t', index_col=0)


def getDNADistMatrix(mat=None, mask_dist=0, gap_dist=0):
    """
    Generates a DNA distance matrix

    Arguments:
    mat = input distance matrix to extend to full alphabet;
          if unspecified, creates Hamming distance matrix that incorporates
          IUPAC equivalencies
    mask_dist = distance for all matches against an N character
    gap_dist = distance for all matches against a gap (-, .) character

    Returns:
    a pandas.DataFrame of distances
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
    mat = input distance matrix to extend to full alphabet;
          if unspecified, creates Hamming distance matrix that incorporates
          IUPAC equivalencies
    mask_dict = score for all matches against an X character
    gap_dist = score for all matches against a gap (-, .) character

    Returns:
    a pandas.DataFrame of distances
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

    :param sequences: list of sequences to be broken into n-mers
    :param n: length of n-mers to return
    :return: dictionary of {sequence: [n-mers]}
    """

    sequences = ['N' * (n-1)/2 + seq + 'N' * (n-1)/2 for seq in sequences]
    nmers = {(seq,[seq[i:i+n] for i in range(len(seq)-n+1)]) for seq in sequences}

    return nmers


def calcDistances(sequences, n, dist_mat, norm):
    """
    Calculate pairwise distances between input sequences

    :param sequences: list of sequences for which to calculate pairwise distances
    :param n: length of n-mers to be used in calculating distance
    :param dist_mat: pandas DataFrame of mutation distances
    :param norm: normalization method
    :return: numpy matrix of pairwise distances between input sequences
    """
    # Initialize output distance matrix
    dists = np.zeros((len(sequences),len(sequences)))
    # Generate dictionary of n-mers from input sequences
    nmers = getNmers(sequences, n)
    # Iterate over combinations of input sequences
    for j,k in combinations(range(len(sequences)), 2):
        # Only consider characters and n-mers with mutations
        mutated = [i for i,(c1,c2) in enumerate(izip(sequences[j],sequences[k])) if c1 != c2]
        seq1 = [sequences[j][i] for i in mutated]
        seq2 = [sequences[k][i] for i in mutated]
        nmer1 = [nmers[seq1][i] for i in mutated]
        nmer2 = [nmers[seq2][i] for i in mutated]

        # Determine normalizing factor
        if norm == 'len':
            norm_by = len(sequences[0])
        elif norm == 'mut':
            norm_by = len(mutated)
        else:
            norm_by = 1

        # Calculate distances
        dists[j,k] = dists[k,j] = \
                sum([dist_mat.at[n2,c1] + dist_mat.at[n1,c2] \
                     for c1,c2,n1,n2 in izip(seq1,seq2,nmer1,nmer2)]) / \
                (2*norm_by)

    return dists


def formClusters(dists, link, distance):
    """
    Form clusters based on hierarchical clustering of input distance matrix with linkage type and cutoff distance
    :param dists: numpy matrix of distances
    :param link: linkage type for hierarchical clustering
    :param distance: distance at which to cut into clusters
    :return: list of cluster assignments
    """
    # Make distance matrix square
    dists = squareform(dists)
    # Compute linkage
    links = linkage(dists, link)
    # Break into clusters based on cutoff
    clusters = fcluster(links, distance, criterion='distance')
    return clusters