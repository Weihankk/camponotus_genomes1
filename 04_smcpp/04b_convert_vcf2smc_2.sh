#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N prep_smc
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=12G
#$ -t 1-31

source activate smcpp_mod

scaffolds=$( head -n${SGE_TASK_ID} vcf_list.txt | tail -n1 )
scaffold_name=${scaffolds%.g.vcf.gz}

smc++ vcf2smc -d 9 9 $scaffolds ${SGE_TASK_ID}_b_modoc.smc $scaffold_name \
modoc:9,12 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 5 5 $scaffolds ${SGE_TASK_ID}_b_laevigatus.smc $scaffold_name \
laevigatus:4,5,10 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 15 15 $scaffolds ${SGE_TASK_ID}_b_herculeanus.smc $scaffold_name \
herculeanus:13,15,16 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 2 2 $scaffolds ${SGE_TASK_ID}_b_pennunk.smc $scaffold_name \
pennunk:1,2,3 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 6 6 $scaffolds ${SGE_TASK_ID}_b_vicinusunk.smc $scaffold_name \
vicinusunk:6,8 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 11 11 $scaffolds ${SGE_TASK_ID}_b_vicinus.smc $scaffold_name \
vicinus:7,11,14 -m ant_smcpp_mask.txt.gz
