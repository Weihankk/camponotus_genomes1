#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N minys1
#$ -q omni
#$ -pe sm 12
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1:17

source activate minys

MinYS.py -1 /lustre/scratch/jmanthey/23_camp1/00_fastq/${SGE_TASK_ID}_R1.fastq.gz \
-2 /lustre/scratch/jmanthey/23_camp1/00_fastq/${SGE_TASK_ID}_R2.fastq.gz \
-ref /lustre/scratch/jmanthey/23_camp1/bloch_penn.fasta \
-out /lustre/scratch/jmanthey/23_camp1/04_minys/${SGE_TASK_ID}
-nb-cores 12
