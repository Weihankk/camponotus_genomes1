#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N smc_ant
#$ -q omni
#$ -pe sm 24
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-6

source activate smcpp_mod

species=$( head -n${SGE_TASK_ID} names_list.txt | tail -n1 )

# run smc++ per species
# 1.983877e-09 mutations / site / year
# generation assumption = 3-4 years for reproductive maturity * 2 = 7.5
# 1.487908e-08 mutations / site / generation

smc++ cv 1.487908e-08 -o ${species}_1 --cores 24 --knots 16 --timepoints 5e3 2e6 --regularization-penalty 6 \
--spline cubic *${species}.smc

# plot output

cd ${species}_1

smc++ plot ${species}.pdf model.final.json -c -g 2
