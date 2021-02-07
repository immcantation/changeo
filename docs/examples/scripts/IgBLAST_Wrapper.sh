#!/usr/bin/env bash
# Download and extract IgBLAST
VERSION="1.17.0"
wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/${VERSION}/ncbi-igblast-${VERSION}-x64-linux.tar.gz
tar -zxf ncbi-igblast-${VERSION}-x64-linux.tar.gz
cp ncbi-igblast-${VERSION}/bin/* ~/bin
# Download reference databases and setup IGDATA directory
fetch_igblastdb.sh -o ~/share/igblast
cp -r ncbi-igblast-${VERSION}/internal_data ~/share/igblast
cp -r ncbi-igblast-${VERSION}/optional_file ~/share/igblast
# Build IgBLAST database from IMGT reference sequences
fetch_imgtdb.sh -o ~/share/germlines/imgt
imgt2igblast.sh -i ~/share/germlines/imgt -o ~/share/igblast
