#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N camp_scaf3
#$ -q omni
#$ -pe sm 24
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

module load gnu

cd working3

~/3d-dna/run-asm-pipeline.sh -r 2 -i 25000 --sort-output \
--editor-coarse-resolution 100000 --editor-coarse-region 300000 --editor-fine-resolution 25000 \
--editor-repeat-coverage 5 \
--polisher-coarse-resolution 100000 --polisher-coarse-region 300000 --polisher-fine-resolution 25000 \
--splitter-coarse-resolution 100000 --splitter-coarse-region 300000 --splitter-fine-resolution 25000 \
~/juicer/references/camponotus_reordered.fasta /lustre/scratch/jmanthey/13_camponotos_genome/scaffolding/aligned/merged_nodups.txt
