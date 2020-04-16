# installed maker already for certhiagenomes1 project

# move to working directory
cd /lustre/scratch/jmanthey/25_camp_maker

# make the control files 
maker -CTL

# edit the control files as necessary
# here, I used four protein files for the first round of maker:

# run maker round 1 script:


#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N camp_r1_maker
#$ -q omni
#$ -pe mpi 108
#$ -P quanah
#$ -l h_rt=48:00:00

mpiexec -n 108 maker -base Camp_r1 maker_r1_opts.ctl maker_bopts.ctl maker_exe.ctl
