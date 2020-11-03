# run repeatmodeler on the camponotus hi-c genome (including all small scaffolds)
# run with a qlogin with 12 processors

# build database of new reference for repeatmodeler
cd /lustre/work/jmanthey/
mkdir camponotus_genome
cd /lustre/work/jmanthey/camponotus_genome

# build the repeatmodeler database
BuildDatabase -engine ncbi -name camp_sp camp_sp_genome_filtered.fasta

# run repeatmodeler
RepeatModeler -database /lustre/work/jmanthey/camponotus_genome/camp_sp -engine ncbi -pa 11



# repeatmodeler for downloaded genome
# build database of new reference for repeatmodeler
cd /lustre/work/jmanthey/camponotus_genome

# build the repeatmodeler database
BuildDatabase -engine ncbi -name form_sely For_selysi_genome.fasta

# run repeatmodeler
RepeatModeler -database /lustre/work/jmanthey/camponotus_genome/form_sely -engine ncbi -pa 11
