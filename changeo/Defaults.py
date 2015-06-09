"""
Default parameters
"""
# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'
from changeo import __version__, __date__

# Imports
import re

# Paths
#default_repo = 'germlines'

# Ig and TCR Regular expressions
allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+[-/\w]*[-\*][\.\w]+))')
gene_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+[-/\w]*))')
family_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+))')

v_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])V[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')
d_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])D[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')
j_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])J[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')

#allele_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*[-\*][\.\w]+)')
#gene_regex = re.compile(r'(IG[HLK][VDJ]\d+[-/\w]*)')
#family_regex = re.compile(r'(IG[HLK][VDJ]\d+)')