"""
Default parameters
"""

# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'

# System settings
default_csv_size = 2**24

# Fields
default_v_field = 'V_CALL'
default_d_field = 'D_CALL'
default_j_field = 'J_CALL'
default_id_field = 'SEQUENCE_ID'
default_seq_field = 'SEQUENCE_IMGT'
default_germ_field = 'GERMLINE_IMGT_D_MASK'
default_junction_field = 'JUNCTION'
default_clone_field = 'CLONE'

# External applications
default_igblast_exec = 'igblastn'
default_tbl2asn_exec = 'tbl2asn'
default_igphyml_exec = 'igphyml'

# Commandline arguments
choices_format = ('changeo', 'airr')
default_format = 'changeo'
default_out_args = {'log_file': None,
                    'out_dir': None,
                    'out_name': None,
                    'out_type': 'tab',
                    'failed': False}
