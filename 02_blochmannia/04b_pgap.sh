#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N bloch_align
#$ -q omni
#$ -pe sm 8
#$ -P quanah
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -t 1-17

source activate ncbi

module load singularity

name_array=$( head -n${SGE_TASK_ID} file_list.txt | tail -n1 )

basename_array=output_${name_array%.fasta}

/lustre/work/jmanthey/ncbi/pgap.py --docker singularity \
-o /lustre/scratch/jmanthey/23_camp1/05_annotate/${basename_array} \
/lustre/scratch/jmanthey/23_camp1/05_annotate/yaml/generic${SGE_TASK_ID}.yaml
