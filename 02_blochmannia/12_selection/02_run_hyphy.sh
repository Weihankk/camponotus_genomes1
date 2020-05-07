#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N hyphy
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1:508

input_array=$( head -n${SGE_TASK_ID} helper1.txt | tail -n1 )

source activate hyphy

hyphy absrel --alignment ${input_array}_aligned_trimmed2.fasta --tree RAxML_bipartitions.${input_array}.tre \
--output absrel___${input_array}.json

hyphy relax --alignment ${input_array}_aligned_trimmed2.fasta \
--tree RAxML_bipartitions.${input_array}.tre.newick.herculeanus --test herculeanus \
--reference background --output relax_herculeanus___${input_array}.json

hyphy relax --alignment ${input_array}_aligned_trimmed2.fasta \
--tree RAxML_bipartitions.${input_array}.tre.newick.laevigatus --test laevigatus \
--reference background --output relax_laevigatus___${input_array}.json

hyphy relax --alignment ${input_array}_aligned_trimmed2.fasta \
--tree RAxML_bipartitions.${input_array}.tre.newick.modoc --test modoc \
--reference background --output relax_modoc___${input_array}.json

hyphy relax --alignment ${input_array}_aligned_trimmed2.fasta \
--tree RAxML_bipartitions.${input_array}.tre.newick.pennunk --test pennunk \
--reference background --output relax_pennunk___${input_array}.json

hyphy relax --alignment ${input_array}_aligned_trimmed2.fasta \
--tree RAxML_bipartitions.${input_array}.tre.newick.vicinus --test vicinus \
--reference background --output relax_vicinus___${input_array}.json

hyphy relax --alignment ${input_array}_aligned_trimmed2.fasta \
--tree RAxML_bipartitions.${input_array}.tre.newick.vicinusunk --test vicinusunk \
--reference background --output relax_vicinusunk___${input_array}.json

