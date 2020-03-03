#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N bloch_align
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-17

module load intel java bwa samtools

/lustre/work/jmanthey/bbmap/bbsplit.sh \
in1=/lustre/scratch/jmanthey/23_camp1/00_fastq/${SGE_TASK_ID}_R1.fastq.gz \
in2=/lustre/scratch/jmanthey/23_camp1/00_fastq/${SGE_TASK_ID}_R2.fastq.gz \
ref=/lustre/scratch/jmanthey/23_camp1/blochmannia_total.fasta \
basename=/lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_%.fastq.gz \
outu1=/lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_R1.fastq.gz \
outu2=/lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_R2.fastq.gz

rm /lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_R1.fastq.gz

rm /lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_R2.fastq.gz

/lustre/work/jmanthey/bbmap/reformat.sh \
in=/lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_blochmannia_total.fastq.gz \
out=/lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_bloch_R1.fastq.gz \
out2=/lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_bloch_R2.fastq.gz

bwa mem -t 2 /lustre/scratch/jmanthey/23_camp1/BX248583_floridanus.fasta \
/lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_bloch_R1.fastq.gz \
/lustre/scratch/jmanthey/23_camp1/05_split/${SGE_TASK_ID}_bloch_R2.fastq.gz > \
/lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}.sam

samtools view -b -S -o /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}.bam \
/lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}.sam

rm /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}.sam

/lustre/work/jmanthey/gatk-4.1.0.0/gatk CleanSam \
-I /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}.bam \
-O /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_cleaned.bam

rm /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}.bam

/lustre/work/jmanthey/gatk-4.1.0.0/gatk SortSam \
-I /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_cleaned.bam \
-O /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_cleaned_sorted.bam --SORT_ORDER coordinate

rm /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_cleaned.bam

/lustre/work/jmanthey/gatk-4.1.0.0/gatk AddOrReplaceReadGroups \
-I /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_cleaned_sorted.bam \
-O /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_cleaned_sorted_rg.bam \
--RGLB 1 --RGPL illumina --RGPU unit1 --RGSM ${SGE_TASK_ID}

rm /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_cleaned_sorted.bam

/lustre/work/jmanthey/gatk-4.1.0.0/gatk MarkDuplicates --REMOVE_DUPLICATES true \
--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 100 \
-M /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_markdups_metric_file.txt \
-I /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_cleaned_sorted_rg.bam \
-O /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_final.bam

rm /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_cleaned_sorted_rg.bam

samtools index /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_final.bam

/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx16g" HaplotypeCaller \
-R /lustre/scratch/jmanthey/23_camp1/BX248583_floridanus.fasta \
-I /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}_final.bam \
-ERC GVCF -O /lustre/scratch/jmanthey/23_camp1/06_bloch_bam/${SGE_TASK_ID}.g.vcf \
--QUIET -ploidy 1

