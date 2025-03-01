blastn -query dna-sequences.fasta -perc_identity 97 -evalue 1e-5 -qcov_hsp_perc 80 -db nt -remote -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" -max_target_seqs 10 -out blastn_results.out

