#!/usr/bin/env python
"""
Presto setup
"""
# Future
from __future__ import division, absolute_import, print_function

# Imports
import os
import sys

# Check setup requirements
if sys.version_info < (2,7,5):
    sys.exit('At least Python 2.7.5 is required.\n')

try:
    from setuptools import setup
except ImportError:
    sys.exit('Please install setuptools before installing changeo.\n')

try:
    from pip.req import parse_requirements
except ImportError:
    sys.exit('Please install pip before installing changeo.\n')

# Get version, author and license information
info_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                         'changeo',
                         'Version.py')
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

# TODO: check pip version to avoid problem with parse_requirements(session=False)
# Parse requirements
try:
    requirements = parse_requirements("requirements.txt", session=False)
except TypeError:
    requirements = parse_requirements("requirements.txt")
install_requires = [str(r.req) for r in requirements]

# Define installation path for commandline tools
scripts = ['AnalyzeAa.py',
           'CreateGermlines.py',
           'DefineClones.py',
           'GapRecords.py',
           'MakeDb.py',
           'ParseDb.py']
install_scripts = [os.path.join('bin', s) for s in scripts]

# Load long package description
with open('README.md', 'r') as f:
    long_description = ''.join([x for x in f])

# Setup
setup(name='changeo',
      version=__version__,
      author=__author__,
      author_email='namita.gupta@yale.edu',
      description='A bioinformatics toolkit for processing high-throughput lymphocyte receptor sequencing data.',
      long_description=long_description,
      zip_safe=False,
      license=__license__,
      url='https://clip.med.yale.edu/changeo',
      keywords='bioinformatics immunoglobulin lymphocyte sequencing',
      install_requires=install_requires,
      packages=['changeo'],
      package_dir={'changeo': 'changeo'},
      package_data={'changeo': ['data/*.tab']},
      scripts=install_scripts,
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'])