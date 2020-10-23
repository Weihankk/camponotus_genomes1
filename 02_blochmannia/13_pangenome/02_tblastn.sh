#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N tblastn
#$ -q omni
#$ -pe sm 12
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

makeblastdb -in blochmannia_sequences_all.fasta -dbtype  nucl -out bloch_all

# output fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
tblastn -query blochmania_AA_all.fasta -out blochmania_AA_all.blast -task tblastn -db bloch_all \
-outfmt 6 -max_target_seqs 100 -evalue 0.1 -num_threads 12
