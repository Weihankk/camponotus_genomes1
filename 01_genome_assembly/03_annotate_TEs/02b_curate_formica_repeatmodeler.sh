# build database of invertebrate repbase + camp laevigatus repeats
cd /lustre/work/jmanthey/
BuildDatabase -engine ncbi -name invert_repbase_test2 invert_repbase_24.03_2.fasta

# rmblast the repeatmodeler output to the new repbase database
# output fields: qseqid qlen sseqid pident length evalue
blastn -query camponotus_genome/form_sely-families.fa -out form_sely_families_repbase.blast -task rmblastn -db invert_repbase_test2 \
-evalue 1e-05 -outfmt '6 qseqid qlen sseqid pident length evalue' -perc_identity 98 -max_target_seqs 10

# remove any repeatmodeler sequences that were >= 98% identical to repbase sequences and > 50% length

# make a blast database of the camponotus genome
cd /lustre/work/jmanthey/camponotus_genome
makeblastdb -in For_selysi_genome.fasta -dbtype nucl -out for_selysi_db1

# blast repeatmodeler output to camponotus genome and take top 100 hits
# output fields: qseqid qlen sseqid sstart send pident length evalue
blastn -query form_sely-families.fasta -out form_sely_repeatmodeler_genomic.blast -task rmblastn -db for_selysi_db1 \
-evalue 1e-20 -outfmt '6 qseqid qlen sseqid sstart send pident length evalue' -perc_identity 90 -qcov_hsp_perc 20 \
-max_target_seqs 100

# filter the top 50 hits for each sequence (if > 50 exist) and create a bed file for extraction of genomic sequence and flanks
# use the r script: bed_creator_from_blast_form.r


# use bedtools to extract fasta sequences from the reference genome based on the previous steps
cd /lustre/work/jmanthey/camponotus_genome/form_bed_extract/
for i in $( ls );
do j=$i.fasta;
bedtools getfasta -fi /lustre/work/jmanthey/camponotus_genome/For_selysi_genome.fasta -bed $i > $j;
done

# take all the extracted fasta files for input into geneious to make mafft alignments and develop consensus sequences

