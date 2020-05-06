#!/usr/bin/env bash
# V segment database
edit_imgt_file.pl IMGT_Human_IGHV.fasta > ~/share/igblast/fasta/imgt_human_ig_v.fasta
makeblastdb -parse_seqids -dbtype nucl -in ~/share/igblast/fasta/imgt_human_ig_v.fasta \
    -out ~/share/igblast/database/imgt_human_ig_v
# D segment database
edit_imgt_file.pl IMGT_Human_IGHD.fasta > ~/share/igblast/fasta/imgt_human_ig_d.fasta
makeblastdb -parse_seqids -dbtype nucl -in ~/share/igblast/fasta/imgt_human_ig_d.fasta \
    -out ~/share/igblast/database/imgt_human_ig_d
# J segment database
edit_imgt_file.pl IMGT_Human_IGHJ.fasta > ~/share/igblast/fasta/imgt_human_ig_j.fasta
makeblastdb -parse_seqids -dbtype nucl -in ~/share/igblast/fasta/imgt_human_ig_j.fasta \
    -out ~/share/igblast/database/imgt_human_ig_j
# Run IgBLAST
export IGDATA=~/share/igblast
igblastn \
    -germline_db_V ~/share/igblast/database/imgt_human_ig_v\
    -germline_db_D ~/share/igblast/database/imgt_human_ig_d \
    -germline_db_J ~/share/igblast/database/imgt_human_ig_v \
    -auxiliary_data ~/share/igblast/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype Ig -organism human \
    -outfmt '7 std qseq sseq btop' \
    -query HD13M.fasta \
    -out HD13M.fmt7
