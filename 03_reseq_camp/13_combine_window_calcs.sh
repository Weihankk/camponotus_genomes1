cd /lustre/scratch/jmanthey/23_camp1/03_vcf/windows

# combine the output for different analyses into a single file each
# first add a header for each file
grep 'pop1' scaffold0001:10000001-10050000__stats.txt > ../window_heterozygosity.txt
grep 'pop1' scaffold0001:10000001-10050000__stats.txt > ../window_titv.txt

# add the relevant stats to each file
for i in $( ls *txt ); do grep 'heterozygosity' $i >> ../window_heterozygosity.txt; done
for i in $( ls *txt ); do grep 'titv' $i >> ../window_titv.txt; done
