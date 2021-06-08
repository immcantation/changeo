"""
Default parameters
"""

# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'

# System settings
default_csv_size = 2**24

# Fields
default_v_field = 'v_call'
default_d_field = 'd_call'
default_j_field = 'j_call'
default_id_field = 'sequence_id'
default_seq_field = 'sequence_alignment'
default_germ_field = 'germline_alignment'
default_junction_field = 'junction'
default_clone_field = 'clone_id'

# Receptor attributes
v_attr = 'v_call'
d_attr = 'd_call'
j_attr = 'j_call'
id_attr = 'sequence_id'
seq_attr = 'sequence_imgt'
germ_attr = 'germline_imgt'
junction_attr = 'junction'
clone_attr = 'clone'

# External applications
default_igblastn_exec = 'igblastn'
default_igblastp_exec = 'igblastp'
default_tbl2asn_exec = 'tbl2asn'
default_igphyml_exec = 'igphyml'

# Commandline arguments
choices_format = ('airr', 'changeo')
default_format = 'airr'
default_out_args = {'log_file': None,
                    'out_dir': None,
                    'out_name': None,
                    'out_type': 'tsv',
                    'failed': False}

# IMGT
default_imgt_id_len = 49
