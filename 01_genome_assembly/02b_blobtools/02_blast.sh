#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N blast_array
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-3707

source activate blast

range_array=$( head -n${SGE_TASK_ID} camp_windows.txt | tail -n1 )

samtools faidx camp_sp_genome.fasta ${range_array} > ${SGE_TASK_ID}.fasta

blastn -query ${SGE_TASK_ID}.fasta -out ${SGE_TASK_ID}.out -task blastn -db nt \
-evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -num_threads 2

rm ${SGE_TASK_ID}.fasta
