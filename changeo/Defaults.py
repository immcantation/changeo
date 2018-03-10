"""
Default parameters
"""

# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'
from changeo import __version__, __date__

# System settings
default_csv_size = 2**24

# Fields
default_v_field = 'V_CALL'
default_d_field = 'D_CALL'
default_j_field = 'J_CALL'
default_clone_field = 'CLONE'

# Commandline arguments
choices_format = ('changeo', 'airr')
default_format = 'changeo'
default_out_args = {'log_file': None,
                    'out_dir': None,
                    'out_name': None,
                    'out_type': None,
                    'failed': True}
