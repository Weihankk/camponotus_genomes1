# some parts of code from https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2


cd /lustre/scratch/jmanthey/25_camp_maker/Camp_r1.maker.output

# get gff from maker without fasta seqs
gff3_merge -n -s -d Camp_r1_master_datastore_index.log > Camp_rnd1.all.maker.noseq.gff

cd ..
mkdir -p augustus/round1
cd augustus/round1

# extract the fasta sequences and surrounding 1000 bp
# some will be missed because bedtools will error when the 1000 bp goes over the edge of the end of the chromosome
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' \
../../Camp_r1.maker.output/Camp_rnd1.all.maker.noseq.gff | \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
bedtools getfasta -fi /home/jmanthey/denovo_genomes/camp_sp_genome_filtered.fasta \
-bed - -fo Camp_rnd1.all.maker.transcripts1000.fasta


# always activate before running (other conda environments have versions of busco that are broken)
source activate busco

# before running busco edit and export the path to the config file
export BUSCO_CONFIG_FILE="/home/jmanthey/busco/config_camp.ini"

run_busco --in Camp_rnd1.all.maker.transcripts1000.fasta \
--out busco_output_round1 --lineage_path /home/jmanthey/busco/hymenoptera_odb9/ \
--mode genome -sp human -c 36 -z --long --augustus_parameters='--progress=true'


# when the busco run finishes (this took 2.5 days on 36 processors)
cd /lustre/scratch/jmanthey/25_camp_maker/augustus/round1/run_busco_output_round1/augustus_output/retraining_parameters
rename 'BUSCO_busco_output_round1_2916219093' 'Camp' *
sed -i 's/BUSCO_busco_output_round1_2916219093/Camp/g' Camp_parameters.cfg
sed -i 's/BUSCO_busco_output_round1_2916219093/Camp/g' Camp_parameters.cfg.orig1 
mkdir $AUGUSTUS_CONFIG_PATH/species/Camp
cp Camp_* $AUGUSTUS_CONFIG_PATH/species/Camp/


