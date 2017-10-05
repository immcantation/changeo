#!/usr/bin/env bash
# Download reference databases
fetch_igblastdb.sh -o ~/share/igblast
fetch_imgtdb.sh -o ~/share/germlines/imgt
# Build IgBLAST database from IMGT reference sequences
imgt2igblast.sh -i ~/share/germlines/imgt -o ~/share/igblast
# Run IgBLAST for human B cell data
run_igblast.sh -s reads.fasta -g human -t ig -b ~/share/igblast