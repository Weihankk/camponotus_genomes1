#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N get
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

update_blastdb.pl nt --decompress
