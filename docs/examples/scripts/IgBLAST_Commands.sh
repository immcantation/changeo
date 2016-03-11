edit_imgt_file.pl IMGT_Human_IGHV.fasta > database/human_igh_v
makeblastdb -parse_seqids -dbtype nucl -in database/human_igh_v

edit_imgt_file.pl IMGT_Human_IGHD.fasta > database/human_igh_d
makeblastdb -parse_seqids -dbtype nucl -in database/human_igh_d

edit_imgt_file.pl IMGT_Human_IGHJ.fasta > database/human_igh_v
makeblastdb -parse_seqids -dbtype nucl -in database/human_igh_j

igblastn \
-germline_db_V database/human_igh_v \
-germline_db_V database/human_igh_v \
-germline_db_V database/human_igh_v \
-auxiliary_data optional_file/human_gl.aux \
-domain_system imgt -ig_seqtype Ig -organism human \
-outfmt '7 std qseq sseq btop' \
-query S43_atleast-2.fasta \
-out S43_atleast-2.fmt7

MakeDb.py igblast -s S43_atleast-2.fasta -i S43_atleast-2.fmt7 -r Human_IGH[VDJ].fasta