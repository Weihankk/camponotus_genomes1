# some of this code modified from https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2

# prepare gff files for round 2 of maker

cd /lustre/scratch/jmanthey/25_camp_maker/Camp_r1.maker.output

gff3_merge -n -s -d Camp_r1_master_datastore_index.log > Camp_round1.all.maker.noseq.gff

# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' Camp_round1.all.maker.noseq.gff > Camp_round1.all.maker.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' Camp_round1.all.maker.noseq.gff > Camp_round1.all.maker.repeats.gff

# move back up to parent directory
cd ..

# copy the original control file with options to save that information and modify the original control file
# for round 2
cp maker_r1_opts.ctl maker_r2_opts.ctl

# rename the maker round 1 output directory
mv Camp_r1.maker.output/ maker_output_round1/

# modify the control file (copied in this directory = maker_r2_opts.ctl)

# run maker round 2 (this took ~ 12 hours)
mpiexec -n 108 maker maker_r2_opts.ctl maker_bopts.ctl maker_exe.ctl

# summarize maker output

cd /lustre/scratch/jmanthey/25_camp_maker/camp_sp_genome_filtered.maker.output

gff3_merge -s -d camp_sp_genome_filtered_master_datastore_index.log > camp_round2.all.maker.gff
fasta_merge -d camp_sp_genome_filtered_master_datastore_index.log
gff3_merge -s  -n -d camp_sp_genome_filtered_master_datastore_index.log > camp_round2.all.maker.noseqs.gff
