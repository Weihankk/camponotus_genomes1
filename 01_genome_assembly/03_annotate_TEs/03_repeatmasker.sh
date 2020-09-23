#Run repeatmasker to annotate where repeats are in the camponotus genome
cd denovo_genomes
RepeatMasker -pa 24 -s -lib ~/RepeatMasker/Libraries/invert_repbase_24.03_custom.fa camp_sp_genome_filtered.fasta
