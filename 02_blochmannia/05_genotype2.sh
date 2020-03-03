#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N bloch_align
#$ -q omni
#$ -pe sm 8
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

module load intel java bwa samtools

/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx20g" GenomicsDBImport \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/1.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/2.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/3.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/4.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/5.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/6.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/7.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/8.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/9.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/10.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/11.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/12.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/13.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/14.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/15.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/16.g.vcf \
-V /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/17.g.vcf \
--genomicsdb-workspace-path /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/bloch_db

/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx64g" GenotypeGVCFs \
-R /lustre/scratch/jmanthey/23_camp1/BX248583_floridanus.fasta \
-V gendb:///lustre/scratch/jmanthey/23_camp1/06_bloch_bam/bloch_db \
--include-non-variant-sites -O /lustre/scratch/jmanthey/23_camp1/07_bloch_vcf/total_bloch.g.vcf
