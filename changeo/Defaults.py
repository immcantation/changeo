"""
Default parameters
"""

# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'
from changeo import __version__, __date__

# System settings
default_csv_size = 2**24

# Annotation parameters
default_delimiter = ('|', '=', ',')
default_separator = default_delimiter[2]

# Commandline arguments
default_format = 'changeo'
default_out_args = {'log_file':None,
                    'delimiter':default_delimiter,
                    'separator':default_separator,
                    'out_dir':None,
                    'out_name':None,
                    'out_type':None,
                    'failed':True}
