#Run repeatmasker to annotate where repeats are in the camponotus genome
cd denovo_genomes

#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N camp_rm
#$ -q omni
#$ -pe sm 24
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

RepeatMasker -pa 24 -s -lib ~/RepeatMasker/Libraries/invert_repbase_24.03_custom.fa camp_sp_genome_filtered.fasta
