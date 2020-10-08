# remove all spaces in fasta files and replace with periods
sed -i 's/ /./g' camponotus_cds_renamed.fasta
sed -i 's/ /./g' lasius_cds_renamed.fasta
sed -i 's/ /./g' formica_cds_renamed.fasta
sed -i 's/ /./g' nylanderia_cds_renamed.fasta

# make blast databases for each of the four species' CDS files

makeblastdb -in camponotus_cds_renamed.fasta -dbtype nucl -out camponotus_cds
makeblastdb -in lasius_cds_renamed.fasta -dbtype nucl -out lasius_cds
makeblastdb -in formica_cds_renamed.fasta -dbtype nucl -out formica_cds
makeblastdb -in nylanderia_cds_renamed.fasta -dbtype nucl -out nylanderia_cds

# blast camponotus to each of the others' databases
# output fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

blastn -query camponotus_cds_renamed.fasta -out q_camponotus_s_lasius.blast -task blastn -db lasius_cds \
-evalue 1e-6 -outfmt 6 -max_target_seqs 100 -num_threads 12

blastn -query camponotus_cds_renamed.fasta -out q_camponotus_s_formica.blast -task blastn -db formica_cds \
-evalue 1e-6 -outfmt 6 -max_target_seqs 100 -num_threads 12

blastn -query camponotus_cds_renamed.fasta -out q_camponotus_s_nylanderia.blast -task blastn -db nylanderia_cds \
-evalue 1e-6 -outfmt 6 -max_target_seqs 100 -num_threads 12

# blast other species to camponotus database

blastn -query lasius_cds_renamed.fasta -out q_lasius_s_camponotus.blast -task blastn -db camponotus_cds \
-evalue 1e-6 -outfmt 6 -max_target_seqs 100 -num_threads 12

blastn -query formica_cds_renamed.fasta -out q_formica_s_camponotus.blast -task blastn -db camponotus_cds \
-evalue 1e-6 -outfmt 6 -max_target_seqs 100 -num_threads 12

blastn -query nylanderia_cds_renamed.fasta -out q_nylanderia_s_camponotus.blast -task blastn -db camponotus_cds \
-evalue 1e-6 -outfmt 6 -max_target_seqs 100 -num_threads 12



# run the r script 04_reciprocal_blast_parsing.r 
# this script looks at all of the blast results and finds the genes where camponotus and the other species have reciprocal
# best hits for camponotus genes
# this then makes a summary table of the names of all of these genes and a fasta file for each group
# for future multiple alignment and tests for selection, etc.
