cd /lustre/scratch/jmanthey/23_camp1/06_cds_output

# make blast databases for camponotus floridanus CDS file

makeblastdb -in BX248583_floridanus_cds.fasta -dbtype  nucl -out cflo_cds

# blast each CDS fasta to the cflo_cds database
# output fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
for i in $( ls *fasta );
do j=${i%fasta}blast;
blastn -query $i -out $j -task blastn -db cflo_cds \
-evalue 1e-2 -outfmt 6 -max_target_seqs 100 -num_threads 2;
done





# copy the 06_cds_output working directory to local computer for processing in R
