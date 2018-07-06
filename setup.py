#!/usr/bin/env python3
"""
Presto setup
"""
# Imports
import os
import sys

# Check setup requirements
if sys.version_info < (3,4,0):
    sys.exit('At least Python 3.4.0 is required.\n')

try:
    from setuptools import setup
except ImportError:
    sys.exit('Please install setuptools before installing changeo.\n')

# Get version, author and license information
info_file = os.path.join('changeo', 'Version.py')
__version__, __author__, __license__ = None, None, None
try:
    exec(open(info_file).read())
except:
    sys.exit('Failed to load package information from %s.\n' % info_file)

if __version__ is None:
    sys.exit('Missing version information in %s\n.' % info_file)
if __author__ is None:
    sys.exit('Missing author information in %s\n.' % info_file)
if __license__ is None:
    sys.exit('Missing license information in %s\n.' % info_file)

# Define installation path for commandline tools
scripts = ['AlignRecords.py',
           'AssignGenes.py',
           'BuildTrees.py',
           'ConvertDb.py',
           'CreateGermlines.py',
           'DefineClones.py',
           'MakeDb.py',
           'ParseDb.py']
install_scripts = [os.path.join('bin', s) for s in scripts]

# Load long package description
#desc_files = ['README.rst', 'INSTALL.rst', 'NEWS.rst']
desc_files = ['README.rst']
long_description = '\n\n'.join([open(f, 'r').read() for f in desc_files])

# Parse requirements
if os.environ.get('READTHEDOCS', None) == 'True':
    # Set empty install_requires to get install to work on readthedocs
    install_requires = []
else:
    with open('requirements.txt') as req:
        install_requires = req.read().splitlines()

# Setup
setup(name='changeo',
      version=__version__,
      author=__author__,
      author_email='immcantation@googlegroups.com',
      description='A bioinformatics toolkit for processing high-throughput lymphocyte receptor sequencing data.',
      long_description=long_description,
      zip_safe=False,
      license=__license__,
      url='http://changeo.readthedocs.io',
      download_url='https://bitbucket.org/kleinstein/changeo/downloads',
      keywords=['bioinformatics', 'sequencing', 'immunology', 'adaptive immunity',
                'immunoglobulin', 'AIRR-seq', 'Rep-Seq',
                'B cell repertoire analysis', 'adaptive immune receptor repertoires'],
      install_requires=install_requires,
      packages=['changeo'],
      package_dir={'changeo': 'changeo'},
      package_data={'changeo': ['data/*.tsv']},
      scripts=install_scripts,
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 3.4',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'])
