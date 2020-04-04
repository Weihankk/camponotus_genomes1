cd /lustre/scratch/jmanthey/23_camp1/03_vcf/

mkdir fasta

cd fasta

mv ../windows/*fasta .

mkdir ../trees

ls *fasta -1 > helper_fasta.txt

cp helper_fasta.txt ../trees/
