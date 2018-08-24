"""
Gene annotations
"""

# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import re
from collections import OrderedDict

# Presto and changeo imports
from changeo.Defaults import default_v_field, default_d_field, default_j_field, default_seq_field

# Ig and TCR Regular expressions
allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+[-/\w]*[-\*][\.\w]+))')
gene_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+[-/\w]*))')
family_regex = re.compile(r'((IG[HLK]|TR[ABGD])([VDJ][A-Z0-9]+))')

v_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])V[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')
d_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])D[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')
j_allele_regex = re.compile(r'((IG[HLK]|TR[ABGD])J[A-Z0-9]+[-/\w]*[-\*][\.\w]+)')

allele_number_regex = re.compile(r'(?<=\*)([\.\w]+)')
c_gene_regex = re.compile(r'((IG[HLK]|TR[ABGD])([DMAGEC][P0-9]?[A-Z]?))')


# TODO:  might be cleaner as getAllele(), getGene(), getFamily()
def parseAllele(alleles, regex, action='first'):
    """
    Extract alleles from strings

    Arguments:
      alleles : string with allele calls
      regex : compiled regular expression for allele match
      action : action to perform for multiple alleles;
               one of ('first', 'set', 'list').
    Returns:
      str : String of the allele when action is 'first';
      tuple : Tuple of allele calls for 'set' or 'list' actions.
    """
    try:
        match = [x.group(0) for x in regex.finditer(alleles)]
    except:
        match = None

    if action == 'first':
        return match[0] if match else None
    elif action == 'set':
        return tuple(sorted(set(match))) if match else None
    elif action == 'list':
        return tuple(sorted(match)) if match else None
    else:
        return None


# TODO: this is not generalized for non-IMGT gapped sequences!
def getVGermline(receptor, references, v_field=default_v_field):
    """
    Extract V allele and germline sequence

    Arguments:
      receptor : Receptor object
      references : dictionary of germline sequences
      v_field : field containing the V allele assignment

    Returns:
      tuple : V allele name, V segment germline sequence
    """
    # Extract V allele call
    vgene = receptor.getVAllele(action='first', field=v_field)

    # Build V segment germline sequence
    if vgene is None:
        try:  vlen = int(receptor.v_germ_length_imgt)
        except (TypeError, ValueError):  vlen = 0
        germ_vseq = 'N' * vlen
    elif vgene in references:
        # Define V germline positions
        try:  vstart = int(receptor.v_germ_start_imgt) - 1
        except (TypeError, ValueError):  vstart = 0
        try:  vlen = int(receptor.v_germ_length_imgt)
        except (TypeError, ValueError):  vlen = 0
        # Define V germline sequence
        vseq = references[vgene]
        vpad = vlen - len(vseq[vstart:])
        if vpad < 0: vpad = 0
        germ_vseq = vseq[vstart:(vstart + vlen)] + ('N' * vpad)
    else:
        germ_vseq = None

    return vgene, germ_vseq


def getDGermline(receptor, references, d_field=default_d_field):
    """
    Extract D allele and germline sequence

    Arguments:
      receptor : Receptor object
      references : dictionary of germline sequences
      d_field : field containing the D allele assignment

    Returns:
      tuple : D allele name, D segment germline sequence
    """
    # Extract D allele call
    dgene = receptor.getDAllele(action='first', field=d_field)

    # Build D segment germline sequence
    if dgene is None:
        germ_dseq = ''
    elif dgene in references:
        # Define D germline positions
        try:  dstart = int(receptor.d_germ_start) - 1
        except (TypeError, ValueError):  dstart = 0
        try:  dlen = int(receptor.d_germ_length)
        except (TypeError, ValueError):  dlen = 0
        # Define D germline sequence
        dseq = references[dgene]
        germ_dseq = dseq[dstart:(dstart + dlen)]
    else:
        germ_dseq = None

    return dgene, germ_dseq


def getJGermline(receptor, references, j_field=default_j_field):
    """
    Extract J allele and germline sequence

    Arguments:
      receptor : Receptor object
      references : dictionary of germline sequences
      j_field : field containing the J allele assignment

    Returns:
      tuple : J allele name, J segment germline sequence
    """
    # Extract J allele call
    jgene = receptor.getJAllele(action='first', field=j_field)

    # Build J segment germline sequence
    if jgene is None:
        try:  jlen = int(receptor.j_germ_length)
        except (TypeError, ValueError):  jlen = 0
        germ_jseq = 'N' * jlen
    elif jgene in references:
        jseq = references[jgene]
        # Define J germline positions
        try:  jstart = int(receptor.j_germ_start) - 1
        except (TypeError, ValueError):  jstart = 0
        try:  jlen = int(receptor.j_germ_length)
        except (TypeError, ValueError):  jlen = 0
        # Define J germline sequence
        jpad = jlen - len(jseq[jstart:])
        if jpad < 0: jpad = 0
        germ_jseq = jseq[jstart:(jstart + jlen)] + ('N' * jpad)
    else:
        germ_jseq = None

    return jgene, germ_jseq


def stitchVDJ(receptor, v_seq, d_seq, j_seq):
    """
    Assemble full length germline sequence

    Arguments:
      receptor : Receptor object
      v_seq : V segment sequence as a string
      d_seq : D segment sequence as a string
      j_seq : J segment sequence as a string

    Returns:
      str : full germline sequence
    """
    # Assemble pieces starting with V segment
    sequence = v_seq

    # Add Ns for first N/P region
    try:  np1_len = int(receptor.np1_length)
    except (TypeError, ValueError):  np1_len = 0
    sequence += 'N' * np1_len

    # Add D segment
    sequence += d_seq

    # Add Ns for second N/P region
    try:  np2_len = int(receptor.np2_length)
    except (TypeError, ValueError):  np2_len = 0
    sequence += 'N' * np2_len

    # Add J segment
    sequence += j_seq

    return sequence


def stitchRegions(receptor, v_seq, d_seq, j_seq):
    """
    Assemble full length region encoding

    Arguments:
      receptor : Receptor object
      v_seq : V segment germline sequence as a string
      d_seq : D segment germline sequence as a string
      j_seq : J segment germline sequence as a string

    Returns:
      str : string defining germline regions
    """
    # Set mode for region definitions
    full_junction = True if getattr(receptor, 'n1_length', None) is not None else False

    # Assemble pieces starting with V segment
    regions = 'V' * len(v_seq)

    # NP nucleotide additions after V
    if not full_junction:
        # PNP nucleotide additions after V
        try:  np1_len = int(receptor.np1_length)
        except (TypeError, ValueError):  np1_len = 0
        regions += 'N' * np1_len
    else:
        # P nucleotide additions before N1
        try:  p3v_len = int(receptor.p3v_length)
        except (TypeError, ValueError):  p3v_len = 0
        # N1 nucleotide additions
        try:  n1_len = int(receptor.n1_length)
        except (TypeError, ValueError):  n1_len = 0
        # P nucleotide additions before D
        try:  p5d_len = int(receptor.p5d_length)
        except (TypeError, ValueError):  p5d_len = 0

        # Update regions
        regions += 'P' * p3v_len
        regions += 'N' * n1_len
        regions += 'P' * p5d_len

    # Add D segment
    regions += 'D' * len(d_seq)

    # NP nucleotide additions before J
    if not full_junction:
        # NP nucleotide additions
        try:  np2_len = int(receptor.np2_length)
        except (TypeError, ValueError):  np2_len = 0
        regions += 'N' * np2_len
    else:
        # P nucleotide additions after D
        try: p3d_len = int(receptor.p3d_length)
        except (TypeError, ValueError): p3d_len = 0
        # N2 nucleotide additions
        try:  n2_len = int(receptor.n2_length)
        except (TypeError, ValueError): n2_len = 0
        # P nucleotide additions before J
        try:  p5j_len = int(receptor.p5j_length)
        except (TypeError, ValueError):  p5j_len = 0

        # Update regions
        regions += 'P' * p3d_len
        regions += 'N' * n2_len
        regions += 'P' * p5j_len

    # Add J segment
    regions += 'J' * len(j_seq)

    return regions


# TODO: Should do 'first' method for ambiguous V/J groups. And explicit allele extraction.
def buildGermline(receptor, references, seq_field=default_seq_field, v_field=default_v_field,
                  d_field=default_d_field, j_field=default_j_field):
    """
    Join gapped germline sequences aligned with sample sequences

    Arguments:
      receptor : Receptor object
      references : dictionary of IMGT gapped germline sequences
      seq_field : field in which to look for sequence
      v_field : field in which to look for V call
      d_field : field in which to look for V call
      j_field : field in which to look for V call

    Returns:
      tuple : log dictionary, dictionary of {germline_type: germline_sequence}, dictionary of {segment: gene call}
    """
    # Return objects
    log = OrderedDict()
    germlines = {'full': '', 'dmask': '', 'vonly': '', 'regions': ''}

    # Build V segment germline sequence
    vgene, germ_vseq = getVGermline(receptor, references, v_field=v_field)
    log['V_CALL'] = vgene
    if germ_vseq is None:
        log['ERROR'] = 'Allele %s in not in the provided germline database.' % vgene
        return log, None, None

    # Build D segment germline sequence
    dgene, germ_dseq = getDGermline(receptor, references, d_field=d_field)
    log['D_CALL'] = dgene
    if germ_dseq is None:
        log['ERROR'] = 'Allele %s in not in the provided germline database.' % vgene
        return log, None, None

    # Build J segment germline sequence
    jgene, germ_jseq = getJGermline(receptor, references, j_field=j_field)
    log['J_CALL'] = jgene
    if germ_jseq is None:
        log['ERROR'] = 'Allele %s in not in the provided germline database.' % vgene
        return log, None, None

    # Stitch complete germlines
    germ_seq = stitchVDJ(receptor, germ_vseq, germ_dseq, germ_jseq)
    regions = stitchRegions(receptor, germ_vseq, germ_dseq, germ_jseq)

    # Update log
    log['SEQUENCE'] = receptor.getField(seq_field)
    log['GERMLINE'] = germ_seq
    log['REGIONS'] = regions

    # Check that input and germline sequence match
    if len(receptor.getField(seq_field)) == 0:
        log['ERROR'] = 'Sequence is missing from the %s field' % seq_field
        return log, None, None

    len_check = len(germ_seq) - len(receptor.getField(seq_field))
    if len_check != 0:
        log['ERROR'] = 'Germline sequence differs in length from input sequence by %i nucleotides.' % abs(len_check)
        return log, None, None

    # Define return germlines object
    germ_dmask = germ_seq[:len(germ_vseq)] + \
                 'N' * (len(germ_seq) - len(germ_vseq) - len(germ_jseq)) + \
                 germ_seq[-len(germ_jseq):]
    germlines = {'full': germ_seq, 'dmask': germ_dmask, 'vonly': germ_vseq, 'regions': regions}
    for k, v in germlines.items():  germlines[k] = v.upper()

    # Define return genes object
    genes = {'v': log['V_CALL'],
             'd': log['D_CALL'],
             'j': log['J_CALL']}

    return log, germlines, genes


def buildClonalGermline(receptors, references, seq_field=default_seq_field,
                        v_field=default_v_field, d_field=default_d_field, j_field=default_j_field):
    """
    Determine consensus clone sequence and create germline for clone

    Arguments:
      receptors : list of Receptor objects
      references : dictionary of IMGT gapped germline sequences
      seq_field : field in which to look for sequence
      v_field : field in which to look for V call
      d_field : field in which to look for D call
      j_field : field in which to look for J call

    Returns:
      tuple : log dictionary, dictionary of {germline_type: germline_sequence},
              dictionary of consensus {segment: gene call}
    """
    # Log
    log = OrderedDict()

    # Create dictionaries to count observed V/J calls
    v_dict = OrderedDict()
    j_dict = OrderedDict()

    # Find longest sequence in clone
    max_length = 0
    for rec in receptors:
        v = rec.getVAllele(action='first', field=v_field)
        v_dict[v] = v_dict.get(v, 0) + 1
        j = rec.getJAllele(action='first', field=j_field)
        j_dict[j] = j_dict.get(j, 0) + 1
        seq_len = len(rec.getField(seq_field))
        if seq_len > max_length:
            max_length = seq_len

    # Consensus V and J having most observations
    v_cons = [k for k in list(v_dict.keys()) if v_dict[k] == max(v_dict.values())]
    j_cons = [k for k in list(j_dict.keys()) if j_dict[k] == max(j_dict.values())]

    # Consensus sequence(s) with consensus V/J calls and longest sequence
    cons = [x for x in receptors if x.getVAllele(action='first', field=v_field) in v_cons and \
                                    x.getJAllele(action='first', field=j_field) in j_cons and \
                                    len(x.getField(seq_field)) == max_length]
    # Consensus sequence(s) with consensus V/J calls but not the longest sequence
    if not cons:
        cons = [x for x in receptors if x.getVAllele(action='first', field=v_field) in v_cons and \
                                        x.getJAllele(action='first', field=j_field) in j_cons]

    # Return without germline if no sequence has both consensus V and J call
    if not cons:
        log['V_CALL'] = ','.join(v_cons)
        log['J_CALL'] = ','.join(j_cons)
        log['ERROR'] = 'No sequence found with both consensus V and J calls.'
        return log, None, None

    # Select consensus Receptor, resolving ties by alphabetical ordering of sequence id.
    cons = sorted(cons, key=lambda x: x.sequence_id)[0]

    # Pad end of consensus sequence with gaps to make it the max length
    gap_length = max_length - len(cons.getField(seq_field))
    if gap_length > 0:
        cons.j_germ_length = int(cons.j_germ_length or 0) + gap_length
        cons.setField(seq_field, cons.getField(seq_field) + ('N' * gap_length))

    # Update lengths padded to longest sequence in clone
    for rec in receptors:
        x = max_length - len(rec.getField(seq_field))
        rec.j_germ_length = int(rec.j_germ_length or 0) + x
        rec.setField(seq_field, rec.getField(seq_field) + ('N' * x))

    # Stitch consensus germline
    cons_log, germlines, genes = buildGermline(cons, references, seq_field=seq_field, v_field=v_field,
                                               d_field=d_field, j_field=j_field)

    # Update log
    log['CONSENSUS'] = cons.sequence_id
    log.update(cons_log)

    # Return log
    return log, germlines, genes