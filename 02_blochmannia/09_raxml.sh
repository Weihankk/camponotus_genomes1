#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N c_raxml
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=12:00:00
#$ -l h_vmem=12G
#$ -t 1-508

# working directory = /lustre/scratch/jmanthey/23_camp1/07_unaligned_fasta
# gene_list == ls -1 *aligned_trimmed* > gene_list.fasta

fasta_array=$( head -n${SGE_TASK_ID} gene_list.fasta | tail -n1 )

input_array=${fasta_array}

output_array=${fasta_array%_aligned_trimmed.fasta}

raxmlHPC-PTHREADS-SSE3 -T 1 -f a -x 50 -m GTRGAMMA -p 253 -N 100 \
-s /lustre/scratch/jmanthey/23_camp1/07_unaligned_fasta/${input_array} \
-n ${output_array}.tre \
-w /lustre/scratch/jmanthey/23_camp1/07_unaligned_fasta/

rm /lustre/scratch/jmanthey/23_camp1/07_unaligned_fasta/RAxML_bestTree.${output_array}.tre

rm /lustre/scratch/jmanthey/23_camp1/07_unaligned_fasta/RAxML_bipartitionsBranchLabels.${output_array}.tre

rm /lustre/scratch/jmanthey/23_camp1/07_unaligned_fasta/RAxML_bootstrap.${output_array}.tre

rm /lustre/scratch/jmanthey/23_camp1/07_unaligned_fasta/RAxML_info.${output_array}.tre
