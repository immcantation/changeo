#!/usr/bin/env python3
"""
Convert simple yaml definitions files to TSV
"""

# Imports
import os, sys
import csv
import yaml
import yamlordereddictloader

# Convert ab1 file
with open('receptor.yaml', 'r') as f:
    records = yaml.load(f, Loader=yaml.FullLoader)

with open('receptor.tsv', 'w') as f:
    writer = csv.writer(f, dialect='excel-tab')
    writer.writerow(['Name', 'Type', 'Description', 'Change-O', 'AIRR', 'IMGT'])
    for r in records['Receptor']:
        writer.writerow([r,
                         records['Receptor'][r]['type'],
                         records['Receptor'][r]['description'],
                         records['Receptor'][r]['changeo'],
                         records['Receptor'][r]['airr'],
                         records['Receptor'][r]['imgt']])


