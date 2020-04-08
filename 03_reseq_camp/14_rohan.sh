#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N ant_rohan
#$ -q omni
#$ -pe sm 8
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1:19

source activate rohan

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jmanthey/anaconda2/envs/rohan/lib/

/lustre/work/jmanthey/rohan/src/rohan --tstv 3.904 -t 8 -o ${SGE_TASK_ID} \
/home/jmanthey/denovo_genomes/camp_sp_genome_filtered.fasta \
/lustre/scratch/jmanthey/23_camp1/01_bam_files/${SGE_TASK_ID}_final.bam


