cd juicer/references

# bwa index
bwa index camponotus_reordered.fasta

cd ../restriction_sites/

# use the juicer python script (edited with arima hi-c restriction sites)
../misc/generate_site_positions2.py Arima camponotus ../references/camponotus_reordered.fasta
