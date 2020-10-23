#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N tblastn
#$ -q omni
#$ -pe sm 4
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

export BATCH_SIZE=10

makeblastdb -in blochmannia_sequences_all.fasta -dbtype  nucl -out bloch_all

tblastn -query blochmania_AA_all.fasta -out blochmania_AA_all.blast -task tblastn -db bloch_all \
-outfmt "6 qseqid sseqid pident length evalue" -max_target_seqs 100 -evalue 0.1 -num_threads 4
