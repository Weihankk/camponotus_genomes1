cd /home/jmanthey/denovo_genomes/

samtools faidx camp_sp_genome_filtered.fasta 

bwa index camp_sp_genome_filtered.fasta

cd

java -jar picard.jar CreateSequenceDictionary R=/home/jmanthey/denovo_genomes/camp_sp_genome_filtered.fasta \
O=/home/jmanthey/denovo_genomes/camp_sp_genome_filtered.dict
