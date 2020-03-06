#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N c_align
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=12:00:00
#$ -l h_vmem=12G
#$ -t 1-542

# working directory = /lustre/scratch/jmanthey/23_camp1/07_unaligned_fasta

fasta_array=$( head -n${SGE_TASK_ID} gene_list.fasta | tail -n1 )

name_array=${fasta_array%.fasta}

# translate the sequences in the fasta file
t_coffee -other_pg seq_reformat -in ${name_array}.fasta -action +translate -output fasta_seq > ${name_array}_protein.fasta

# align the sequences
t_coffee ${name_array}_protein.fasta -mode mcoffee

# back translate to nucleotides
t_coffee -other_pg seq_reformat -in ${name_array}.fasta -in2 ${name_array}_protein.aln -action +thread_dna_on_prot_aln \
-output fasta_aln > ${name_array}_aligned.fasta

# remove gaps from alignment
trimal -in ${name_array}_aligned.fasta -out ${name_array}_aligned_trimmed.fasta -nogaps -fasta
