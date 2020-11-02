# run the msmc utility script for each directory to create bootstrap replicates
cd /lustre/scratch/jmanthey/23_camp1/msmc_demography

for i in $( ls ); do 
cd $i;
~/multihetsep_bootstrap.py -n 10 -s 1000000 --chunks_per_chromosome 11 --nr_chromosomes 31 \
--seed 324324 bootstrap *txt;
cd ..;
done
