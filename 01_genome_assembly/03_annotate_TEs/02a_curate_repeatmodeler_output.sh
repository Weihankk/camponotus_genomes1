# build database of vertebrate repbase 
cd /lustre/work/jmanthey/
BuildDatabase -engine ncbi -name invert_repbase invert_repbase_24.03.fasta

# rmblast the repeatmodeler output to the new repbase database
# output fields: qseqid qlen sseqid pident length evalue
blastn -query camponotus_genome/camp_sp-families.fa -out camp_sp_families_repbase.blast -task rmblastn -db invert_repbase \
-evalue 1e-05 -outfmt '6 qseqid qlen sseqid pident length evalue' -perc_identity 98 -max_target_seqs 10

# remove any repeatmodeler sequences that were >= 98% identical to repbase sequences and > 50% length

# make a blast database of the camponotus genome
cd /lustre/work/jmanthey/camponotus_genome
makeblastdb -in camp_sp_genome_filtered.fasta -dbtype nucl -out camp_sp_db1

# blast repeatmodeler output to camponotus genome and take top 100 hits
# output fields: qseqid qlen sseqid sstart send pident length evalue
blastn -query camp_sp-families.fasta -out camp_sp_repeatmodeler_genomic.blast -task rmblastn -db camp_sp_db1 \
-evalue 1e-20 -outfmt '6 qseqid qlen sseqid sstart send pident length evalue' -perc_identity 90 -qcov_hsp_perc 20 \
-max_target_seqs 100

# filter the top 50 hits for each sequence (if > 50 exist) and create a bed file for extraction of genomic sequence and flanks
# use the r script: bed_creator_from_blast.r


# use bedtools to extract fasta sequences from the reference genome based on the previous steps
cd /lustre/work/jmanthey/camponotus_genome/camp_bed_extract/
for i in $( ls );
do j=$i.fasta;
bedtools getfasta -fi /lustre/work/jmanthey/camponotus_genome/camp_sp_genome_filtered.fasta -bed $i > $j;
done

# take all the extracted fasta files for input into geneious to make mafft alignments and develop consensus sequences





# some of the TE consensus sequences needed to be extended further, redo the above steps with the updated fasta
# file that has already been extended. Some are also just for clarification of endpoints of TEs.

# blast repeatmodeler output to certhia genome and take top 100 hits
# output fields: qseqid qlen sseqid sstart send pident length evalue
cd /lustre/work/jmanthey/camponotus_genome
blastn -query refine2_campsp.fasta -out camp_sp_repeatmodeler_genomic2.blast -task rmblastn -db camp_sp_db1 \
-evalue 1e-20 -outfmt '6 qseqid qlen sseqid sstart send pident length evalue' -perc_identity 90 -qcov_hsp_perc 20 \
-max_target_seqs 100

# filter the top 50 hits for each sequence (if > 50 exist) and create a bed file for extraction of genomic sequence and flanks
# use the r script: bed_creator_from_blast2.r

# use bedtools to extract fasta sequences from the reference genome based on the previous steps
cd /lustre/work/jmanthey/camponotus_genome/camp_bed_extract2/
for i in $( ls );
do j=$i.fasta;
bedtools getfasta -fi /lustre/work/jmanthey/camponotus_genome/camp_sp_genome_filtered.fasta -bed $i > $j;
done

# take all the extracted fasta files for input into geneious to make mafft alignments and develop consensus sequences
