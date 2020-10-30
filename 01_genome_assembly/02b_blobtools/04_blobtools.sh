./blobtools create \
 -i camp/camp_sp_genome.fasta \
 -b camp/1_final.bam \
 -t camp/total_blast.out \
 -o camp/camp_blobplot

./blobtools view -r order -i camp/camp_blobplot.blobDB.json -o camp

grep '^##' camp.camp_blobplot.blobDB.table.txt ; \
grep -v '^##' camp.camp_blobplot.blobDB.table.txt | \
column -t -s $'\t'
 
./blobtools plot \
-i camp/camp_blobplot.blobDB.json \
-o camp/
