#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N camp_depth
#$ -q omni
#$ -pe sm 4
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=12G

samtools depth -a 1_final.bam 2_final.bam 3_final.bam 4_final.bam 5_final.bam 6_final.bam 7_final.bam 8_final.bam \
9_final.bam 10_final.bam 11_final.bam 12_final.bam 13_final.bam 14_final.bam 15_final.bam 16_final.bam 17_final.bam \
18_final.bam 19_final.bam > camp_coverage.txt
