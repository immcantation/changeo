#!/usr/bin/env bash
# V-segment database
edit_imgt_file.pl IMGT_Human_IGHV.fasta > database/human_igh_v
makeblastdb -parse_seqids -dbtype nucl -in database/human_igh_v
# D-segment database
edit_imgt_file.pl IMGT_Human_IGHD.fasta > database/human_igh_d
makeblastdb -parse_seqids -dbtype nucl -in database/human_igh_d
# J-segment database
edit_imgt_file.pl IMGT_Human_IGHJ.fasta > database/human_igh_j
makeblastdb -parse_seqids -dbtype nucl -in database/human_igh_j
# Run IgBLAST
igblastn \
    -germline_db_V database/human_igh_v \
    -germline_db_D database/human_igh_d \
    -germline_db_J database/human_igh_j \
    -auxiliary_data optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype Ig -organism human \
    -outfmt '7 std qseq sseq btop' \
    -query S43_atleast-2.fasta \
    -out S43_atleast-2.fmt7
# Parse IgBLAST output
MakeDb.py igblast -i S43_atleast-2.fmt7 -s S43_atleast-2.fasta \
    -r IMGT_Human_IGHV.fasta IMGT_Human_IGHD.fasta IMGT_Human_IGHJ.fasta \
    --regions --scores