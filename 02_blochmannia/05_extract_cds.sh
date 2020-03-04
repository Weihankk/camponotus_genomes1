cd /Volumes/ants_2020/03_blochmannia_analysis/01_extract_cds

rm fasta/*fai
# all whole genomes in fasta directory
for i in $( ls fasta ); do
j=${i%.fasta};
echo $j;
cat output_${j}/annot.gff | grep "CDS" > gff_trimmed/${j}.gff;
bedtools getfasta -s -fi fasta/$i -bed gff_trimmed/${j}.gff > cds_output/${j}.fasta;
done

