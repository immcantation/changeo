"""
Distance data structures and functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'
from changeo import __version__, __date__

# Imports
import pandas as pd
from itertools import product
from pkg_resources import resource_stream

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
    m3n_distance = pd.read_csv(f, sep='\t', index_col=0).to_dict()

with resource_stream(__name__, 'data/HS5F_Distance.tab') as f:
    hs5f_distance = pd.read_csv(f, sep='\t', index_col=0).to_dict()


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
            dist_mat.loc[i, j] = mat.loc[i, j]
    # If no input matrix, create IUPAC-defined Hamming distance
    else:
        for i,j in product(dist_mat.index, dist_mat.columns):
            dist_mat.loc[i, j] = 1 - scoreDNA(i, j,
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
            dist_mat.loc[i, j] = mat.loc[i, j]
    # If no input matrix, create IUPAC-defined Hamming distance
    else:
        for i,j in product(dist_mat.index, dist_mat.columns):
            dist_mat.loc[i, j] = 1 - scoreAA(i, j,
                                            mask_score=(1 - mask_dist, 1 - mask_dist),
                                            gap_score=(1 - gap_dist, 1 - gap_dist))

    return dist_mat