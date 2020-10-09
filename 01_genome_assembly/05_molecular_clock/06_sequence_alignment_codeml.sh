#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N c_align_paml
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=12:00:00
#$ -l h_vmem=12G
#$ -t 1-6940

# working directory = /lustre/scratch/jmanthey/25_camp_maker/cds_fasta_files/camp_sp_genome_filtered.maker.output

# translate the sequences in the fasta file
t_coffee -other_pg seq_reformat -in ${SGE_TASK_ID}.fasta -action +translate -output fasta_seq > ${SGE_TASK_ID}_protein.fasta

# align the sequences
t_coffee ${SGE_TASK_ID}_protein.fasta -mode mcoffee

# back translate to nucleotides
t_coffee -other_pg seq_reformat -in ${SGE_TASK_ID}.fasta -in2 ${SGE_TASK_ID}_protein.aln -action +thread_dna_on_prot_aln \
-output fasta_aln > ${SGE_TASK_ID}_aligned.fasta

# remove gaps from alignment
trimal -in ${SGE_TASK_ID}_aligned.fasta -out ${SGE_TASK_ID}_aligned_trimmed.fasta -nogaps -fasta

# create a new directory for this alignment and move to it
mkdir ${SGE_TASK_ID}_work
cd ${SGE_TASK_ID}_work/

# copy a control file to this directory
cp ../../codeml.ctl .

# copy the alignment to this directory
cp ../${SGE_TASK_ID}_aligned_trimmed.fasta .

# change the control file input and output names 
sed -i "s/blank.fasta/${SGE_TASK_ID}_aligned_trimmed.fasta/g" codeml.ctl
sed -i "s/blank_output.txt/${SGE_TASK_ID}_output.txt/g" codeml.ctl

# run 1st iteration of codeml
codeml codeml.ctl

# copy camponotus control file to this directory
cp ../../codeml_2.ctl .

# change the control file input and output names for camponotus
sed -i "s/blank.fasta/${SGE_TASK_ID}_aligned_trimmed.fasta/g" codeml_2.ctl
sed -i "s/blank_output.txt/${SGE_TASK_ID}_output_camponotus.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree.tre/unrooted_tree_camponotus.tre/g' codeml_2.ctl

# run camponotus iteration of codeml
codeml codeml_2.ctl

# change the control file input and output names for formica
sed -i "s/${SGE_TASK_ID}_output_camponotus.txt/${SGE_TASK_ID}_output_formica.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree_camponotus.tre/unrooted_tree_formica.tre/g' codeml_2.ctl

# run formica iteration of codeml
codeml codeml_2.ctl

# change the control file input and output names for nylanderia
sed -i "s/${SGE_TASK_ID}_output_formica.txt/${SGE_TASK_ID}_output_nylanderia.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree_formica.tre/unrooted_tree_nylanderia.tre/g' codeml_2.ctl

# run nylanderia iteration of codeml
codeml codeml_2.ctl

# change the control file input and output names for lasius
sed -i "s/${SGE_TASK_ID}_output_nylanderia.txt/${SGE_TASK_ID}_output_lasius.txt/g" codeml_2.ctl
sed -i 's/unrooted_tree_nylanderia.tre/unrooted_tree_lasius.tre/g' codeml_2.ctl

# run nylanderia iteration of codeml
codeml codeml_2.ctl

# add all likelihoods and omegas to a single file
grep '^lnL' ${SGE_TASK_ID}_output.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt
grep '^omega' ${SGE_TASK_ID}_output.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_camponotus.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_camponotus.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_formica.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_formica.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_nylanderia.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_nylanderia.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt
grep '^lnL' ${SGE_TASK_ID}_output_lasius.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt
grep '^w (dN/dS)' ${SGE_TASK_ID}_output_lasius.txt >> /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output/output/${SGE_TASK_ID}_total_output.txt

# move the alignment fasta to the output directory 
mv ${SGE_TASK_ID}_aligned_trimmed.fasta ../../output

# remove temp working directory
cd ..
rm -r ${SGE_TASK_ID}_work/
