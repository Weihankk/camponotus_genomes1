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

smc++ vcf2smc -d 12 12 $scaffolds ${SGE_TASK_ID}_modoc.smc $scaffold_name \
modoc:9,12 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 10 10 $scaffolds ${SGE_TASK_ID}_laevigatus.smc $scaffold_name \
laevigatus:4,5,10 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 16 16 $scaffolds ${SGE_TASK_ID}_herculeanus.smc $scaffold_name \
herculeanus:13,15,16 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 3 3 $scaffolds ${SGE_TASK_ID}_pennunk.smc $scaffold_name \
pennunk:1,2,3 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 8 8 $scaffolds ${SGE_TASK_ID}_vicinusunk.smc $scaffold_name \
vicinusunk:6,8 -m ant_smcpp_mask.txt.gz

smc++ vcf2smc -d 14 14 $scaffolds ${SGE_TASK_ID}_vicinus.smc $scaffold_name \
vicinus:7,11,14 -m ant_smcpp_mask.txt.gz
