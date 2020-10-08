cd /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output

# extract the cds for the entire gff
cat camp_round2.all.maker.noseqs.gff | awk '{ if ($3 == "CDS") print $0 }' > camp_round2_cds.gff

# extract sequences (make sure to force strandedness)
bedtools getfasta -s -fi /home/jmanthey/denovo_genomes/camp_sp_genome_filtered.fasta \
-bed camp_round2_cds.gff > camponotus_cds.fasta


# other ant species

# extract the cds for the entire gff
cat /home/jmanthey/references/GCA_001045655.1_ASM104565v1_genomic.gff | awk '{ if ($3 == "CDS") print $0 }' > lasius_cds.gff
cat /home/jmanthey/references/GCF_003651465.1_ASM365146v1_genomic.gff | awk '{ if ($3 == "CDS") print $0 }' > formica_cds.gff
cat /home/jmanthey/references/GCF_005281655.1_TAMU_Nfulva_1.0_genomic.gff | awk '{ if ($3 == "CDS") print $0 }' > nylanderia_cds.gff

# extract sequences (make sure to force strandedness)
bedtools getfasta -s -fi /home/jmanthey/references/GCA_001045655.1_ASM104565v1_genomic.fna \
-bed lasius_cds.gff > lasius_cds.fasta

bedtools getfasta -s -fi /home/jmanthey/references/GCF_003651465.1_ASM365146v1_genomic.fna \
-bed formica_cds.gff > formica_cds.fasta

bedtools getfasta -s -fi /home/jmanthey/references/GCF_005281655.1_TAMU_Nfulva_1.0_genomic.fna \
-bed nylanderia_cds.gff > nylanderia_cds.fasta


# run the r script 02_label_cds.r to rename the fasta file sequences using the gene names in the gff as well as concatenate 
# the genes' multiple cds
