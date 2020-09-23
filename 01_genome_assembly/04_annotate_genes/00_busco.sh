source activate busco

cd /home/jmanthey/denovo_genomes

# before running busco edit and export the path to the config file
export BUSCO_CONFIG_FILE="/home/jmanthey/busco/config.ini"

run_busco --in /home/jmanthey/denovo_genomes/camp_sp_genome_filtered.fasta \
--out  camp_busco --lineage_path /home/jmanthey/busco/hymenoptera_odb9/ \
--mode genome -sp human -c 16 
