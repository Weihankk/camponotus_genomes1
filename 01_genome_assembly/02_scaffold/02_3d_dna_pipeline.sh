# step 1
# run juicer

# Run juicer.sh with the flags -z <path to genome fasta file>, 
# -p <path to mygenome.chrom.sizes>, and -y <path to mygenome_myenzyme.txt>

# need to edit the juicer.sh script to include the project id  '-P quanah'
# remove all 'use' commands from the script and remove all 'load' commands

/home/jmanthey/juicer/scripts/juicer.sh -z /home/jmanthey/juicer/references/camponotus_reordered.fasta \
-p /home/jmanthey/juicer/restriction_sites/camponotus_reordered.contig.sizes \
-y /home/jmanthey/juicer/restriction_sites/camponotus_Arima.txt \
-d /lustre/scratch/jmanthey/13_camponotos_genome/scaffolding \
-q omni -l omni -s Arima -a 'Camponotus genome assembly Hi-C' \
-D juicer -Q 48:00:00 -L 48:00:00 -t 8

# fails after merging alignments and splitting that
# needs to run 02b_merge_no_dupes.r to remove all duplicates from alignments

# combine output from r script for use in the 3d dna pipeline
for i in {1..332}; do 
	cat no_dupes_${i}.txt >> merged_nodups.txt;
done
