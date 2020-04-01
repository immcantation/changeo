"""
Alignment manipulation
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import re
from Bio.Seq import Seq

# Presto and changeo imports
from changeo.Gene import getVAllele, getJAllele

# Load regions
# import yaml
# from pkg_resources import resource_stream
# with resource_stream(__name__, 'data/regions.yaml') as f:
#     imgt_regions = yaml.load(f, Loader=yaml.FullLoader)
imgt_regions = {'default': {'fwr1': 1,
                            'cdr1': 27,
                            'fwr2': 39,
                            'cdr2': 56,
                            'fwr3': 66,
                            'cdr3': 105},
                'rhesus-igl': {'fwr1': 1,
                               'cdr1': 28,
                               'fwr2': 40,
                               'cdr2': 59,
                               'fwr3': 69,
                               'cdr3': 108}}


class RegionDefinition:
    """
    FWR and CDR region boundary definitions
    """
    def __init__(self, junction_length, amino_acid=False, definition='default'):
        """
        Initializer

        Arguments:
          junction_length (int): length of the junction region. If None then CDR3 end and FWR4 start/end are undefined.
          definition (str): region definition entry in the data/regions.yaml file to use.
          amino_acid (bool): if True define boundaries in amino acid space, otherwise use nucleotide positions.

        Returns:
          changeo.Alignment.RegionDefinition
        """
        self.junction_length = junction_length
        self.amino_acid = amino_acid
        self.definition = definition
        pos_mod = 1 if amino_acid else 3

        # Define regions
        regions = {k: (int(v) - 1) * pos_mod for k, v in imgt_regions[definition].items()}

        # Assign positions
        if junction_length is not None:
            fwr4_start = max(regions['cdr3'], regions['cdr3'] - (2 * pos_mod) + junction_length) \
                if junction_length is not None else None
            junction_end = fwr4_start + (1 * pos_mod)
        else:
            fwr4_start = None
            junction_end = None
        self.positions = {'fwr1': [regions['fwr1'], regions['cdr1']],
                          'cdr1': [regions['cdr1'], regions['fwr2']],
                          'fwr2': [regions['fwr2'], regions['cdr2']],
                          'cdr2': [regions['cdr2'], regions['fwr3']],
                          'fwr3': [regions['fwr3'], regions['cdr3']],
                          'cdr3': [regions['cdr3'], fwr4_start],
                          'fwr4': [fwr4_start, None],
                          'junction': [regions['cdr3'] - (1 * pos_mod), junction_end]}

    def getRegions(self, seq):
        """
        Return IMGT defined FWR and CDR regions

        Arguments:
          seq : IMGT-gapped sequence.

        Returns:
          dict : dictionary of FWR and CDR sequences.
        """
        regions = {'fwr1_imgt': None,
                   'fwr2_imgt': None,
                   'fwr3_imgt': None,
                   'fwr4_imgt': None,
                   'cdr1_imgt': None,
                   'cdr2_imgt': None,
                   'cdr3_imgt': None}
        try:
            seq_len = len(seq)
            regions['fwr1_imgt'] = seq[self.positions['fwr1'][0]:min(self.positions['fwr1'][1], seq_len)]
        except (KeyError, IndexError, TypeError):
            return regions

        try:
            regions['cdr1_imgt'] = seq[self.positions['cdr1'][0]:min(self.positions['cdr1'][1], seq_len)]
        except (IndexError):
            return regions

        try:
            regions['fwr2_imgt'] = seq[self.positions['fwr2'][0]:min(self.positions['fwr2'][1], seq_len)]
        except (IndexError):
            return regions

        try:
            regions['cdr2_imgt'] = seq[self.positions['cdr2'][0]:min(self.positions['cdr2'][1], seq_len)]
        except (IndexError):
            return regions

        try:
            regions['fwr3_imgt'] = seq[self.positions['fwr3'][0]:min(self.positions['fwr3'][1], seq_len)]
        except (IndexError):
            return regions

        try:
            regions['cdr3_imgt'] = seq[self.positions['cdr3'][0]:min(self.positions['cdr3'][1], seq_len)]
            regions['fwr4_imgt'] = seq[self.positions['fwr4'][0]:]
        except (KeyError, IndexError, TypeError):
            return regions

        return regions


def decodeBTOP(btop):
    """
    Parse a BTOP string into a list of tuples in CIGAR annotation.

    Arguments:
      btop : BTOP string.

    Returns:
      list : tuples of (operation, length) for each operation in the BTOP string using CIGAR annotation.
    """
    # Determine chunk type and length
    def _recode(m):
        if m.isdigit():  return ('=', int(m))
        elif m[0] == '-':  return ('I', len(m) // 2)
        elif m[1] == '-':  return ('D', len(m) // 2)
        else:  return ('X', len(m) // 2)

    # Split BTOP string into sections
    btop_split = re.sub(r'(\d+|[-A-Z]{2})', r'\1;', btop)
    # Parse each chunk of encoding
    matches = re.finditer(r'(\d+)|([A-Z]{2};)+|(-[A-Z];)+|([A-Z]-;)+', btop_split)

    return [_recode(m.group().replace(';', '')) for m in matches]


def decodeCIGAR(cigar):
    """
    Parse a CIGAR string into a list of tuples.

    Arguments:
      cigar : CIGAR string.

    Returns:
      list : tuples of (operation, length) for each operation in the CIGAR string.
    """
    matches = re.findall(r'(\d+)([A-Z])', cigar)

    return [(m[1], int(m[0])) for m in matches]


def encodeCIGAR(alignment):
    """
    Encodes a list of tuple with alignment information into a CIGAR string.

    Arguments:
      tuple : tuples of (type, length) for each alignment operation.

    Returns:
      str : CIGAR string.
    """
    return ''.join(['%i%s' % (x, s) for s, x in alignment])


def padAlignment(alignment, q_start, r_start):
    """
    Pads the start of an alignment based on query and reference positions.

    Arguments:
      alignment : tuples of (operation, length) for each alignment operation.
      q_start : query (input) start position (0-based)
      r_start : reference (subject) start position (0-based)

    Returns:
      list : updated list of tuples of (operation, length) for the alignment.
    """
    # Copy list to avoid weirdness
    result = alignment[:]

    # Add query deletions
    if result [0][0] == 'S':
        result[0] = ('S', result[0][1] + q_start)
    elif q_start > 0:
        result.insert(0, ('S', q_start))

    # Add reference padding if present
    if result[0][0] == 'N':
        result[0] = ('N', result[0][1] + r_start)
    elif result [0][0] == 'S' and result[1][0] == 'N':
        result[1] = ('N', result[1][1] + r_start)
    elif result[0][0] == 'S' and r_start > 0:
        result.insert(1, ('N', r_start))
    elif r_start > 0:
        result.insert(0, ('N', r_start))

    return result


def alignmentPositions(alignment):
    """
    Extracts start position and length from an alignment

    Arguments:
      alignment : tuples of (operation, length) for each alignment operation.

    Returns:
      dict : query (q) and reference (r) start (0-based) and length information with keys
             {q_start, q_length, r_start, r_length}.
    """
    # Return object
    result = {'q_start': 0,
              'q_length': 0,
              'r_start': 0,
              'r_length': 0}

    # Query start
    if alignment[0][0] == 'S':
        result['q_start'] = alignment[0][1]

    # Reference start
    if alignment[0][0] == 'N':
        result['r_start'] = alignment[0][1]
    elif alignment[0][0] == 'S' and alignment[1][0] == 'N':
        result['r_start'] = alignment[1][1]

    # Reference length
    for x, i in alignment:
        if x in ('M', '=', 'X'):
            result['r_length'] += i
            result['q_length'] += i
        elif x == 'D':
            result['r_length'] += i
        elif x == 'I':
            result['q_length'] += i

    return result


def gapV(seq, v_germ_start, v_germ_length, v_call, references, asis_calls=False):
    """
    Construction IMGT-gapped V segment sequences.

    Arguments:
      seq (str): V(D)J sequence alignment (SEQUENCE_VDJ).
      v_germ_start (int): start position V segment alignment in the germline (V_GERM_START_VDJ, 1-based).
      v_germ_length (int): length of the V segment alignment against the germline (V_GERM_LENGTH_VDJ, 1-based).
      v_call (str): V segment allele assignment (V_CALL).
      references (dict): dictionary of IMGT-gapped reference sequences.
      asis_calls (bool): if True do not parse v_call for allele names and just split by comma.

    Returns:
      dict: dictionary containing IMGT-gapped query sequences and germline positions.

    Raises:
      KeyError: raised if the v_call is not found in the reference dictionary.
    """
    # Initialize return object
    imgt_dict = {'sequence_imgt': None,
                 'v_germ_start_imgt': None,
                 'v_germ_length_imgt': None}

    # Initialize imgt gapped sequence
    seq_imgt = '.' * (int(v_germ_start) - 1) + seq

    # Extract first V call
    if not asis_calls:
        vgene = getVAllele(v_call, action='first')
    else:
        vgene = v_call.split(',')[0]

    # Find gapped germline V segment
    try:
    #if vgene in references:
        vgap = references[vgene]
        # Iterate over gaps in the germline segment
        gaps = re.finditer(r'\.', vgap)
        gapcount = int(v_germ_start) - 1
        for gap in gaps:
            i = gap.start()
            # Break if gap begins after V region
            if i >= v_germ_length + gapcount:
                break
            # Insert gap into IMGT sequence
            seq_imgt = seq_imgt[:i] + '.' + seq_imgt[i:]
            # Update gap counter
            gapcount += 1

        imgt_dict['sequence_imgt'] = seq_imgt
        # Update IMGT positioning information for V
        imgt_dict['v_germ_start_imgt'] = 1
        imgt_dict['v_germ_length_imgt'] = v_germ_length + gapcount
    except KeyError as e:
        raise KeyError('%s was not found in the germline repository.' % vgene)
    #else:
    #    printWarning('%s was not found in the germline repository. IMGT-gapped sequence cannot be determined.' % vgene)

    return imgt_dict


def inferJunction(seq, j_germ_start, j_germ_length, j_call, references, asis_calls=False, regions='default'):
    """
    Identify junction region by IMGT definition.

    Arguments:
      seq (str): IMGT-gapped V(D)J sequence alignment (SEQUENCE_IMGT).
      j_germ_start (int): start position J segment alignment in the germline (J_GERM_START, 1-based).
      j_germ_length (int): length of the J segment alignment against the germline (J_GERM_LENGTH).
      j_call (str): J segment allele assignment (J_CALL).
      references (dict): dictionary of IMGT-gapped reference sequences.
      asis_calls (bool): if True do not parse V_CALL for allele names and just split by comma.
      regions (str): name of the IMGT FWR/CDR region definitions to use.

    Returns:
      dict : dictionary containing junction sequence, translation and length.
    """
    junc_dict = {'junction': None,
                 'junction_aa': None,
                 'junction_length': None}

    # Find germline J segment
    if not asis_calls:
        jgene = getJAllele(j_call, action='first')
    else:
        jgene = j_call.split(',')[0]
    jgerm = references.get(jgene, None)

    if jgerm is not None:
        # Look for (F|W)GXG amino acid motif in germline nucleotide sequence
        motif = re.search(r'T(TT|TC|GG)GG[ACGT]{4}GG[AGCT]', jgerm)

        # Define junction end position
        seq_len = len(seq)
        if motif:
            j_start = seq_len - j_germ_length
            motif_pos = max(motif.start() - j_germ_start + 1, -1)
            junc_end = j_start + motif_pos + 3
        else:
            junc_end = seq_len

        # Extract junction
        rd = RegionDefinition(None, amino_acid=False, definition=regions)
        junc_start = rd.positions['junction'][0]
        junc_dict['junction'] = seq[junc_start:junc_end]
        junc_len = len(junc_dict['junction'])
        junc_dict['junction_length'] = junc_len

        # Translation
        junc_tmp = junc_dict['junction'].replace('-', 'N').replace('.', 'N')
        if junc_len % 3 > 0:  junc_tmp = junc_tmp[:junc_len - junc_len % 3]
        junc_dict['junction_aa'] = str(Seq(junc_tmp).translate())

    return junc_dict


def getRegions(seq, junction_length):
    """
    Identify FWR and CDR regions by IMGT definition.

    Arguments:
      seq : IMGT-gapped sequence.
      junction_length : length of the junction region in nucleotides.

    Returns:
      dict : dictionary of FWR and CDR sequences.
    """
    region_dict = {'fwr1_imgt': None,
                   'fwr2_imgt': None,
                   'fwr3_imgt': None,
                   'fwr4_imgt': None,
                   'cdr1_imgt': None,
                   'cdr2_imgt': None,
                   'cdr3_imgt': None}
    try:
        seq_len = len(seq)
        region_dict['fwr1_imgt'] = seq[0:min(78, seq_len)]
    except (KeyError, IndexError, TypeError):
        return region_dict

    try: region_dict['cdr1_imgt'] = seq[78:min(114, seq_len)]
    except (IndexError): return region_dict

    try: region_dict['fwr2_imgt'] = seq[114:min(165, seq_len)]
    except (IndexError): return region_dict

    try: region_dict['cdr2_imgt'] = seq[165:min(195, seq_len)]
    except (IndexError): return region_dict

    try: region_dict['fwr3_imgt'] = seq[195:min(312, seq_len)]
    except (IndexError): return region_dict

    try:
        # CDR3
        cdr3_end = 306 + junction_length
        region_dict['cdr3_imgt'] = seq[312:cdr3_end]
        # FWR4
        region_dict['fwr4_imgt'] = seq[cdr3_end:]
    except (KeyError, IndexError, TypeError):
        return region_dict

    return region_dict
