# rename the samples for compatability with the cluster array jobs

cd /lustre/scratch/jmanthey/23_camp1/00_fastq/
while read -r name1 name2; do
	mv $name1 $name2
done < rename_samples.txt
