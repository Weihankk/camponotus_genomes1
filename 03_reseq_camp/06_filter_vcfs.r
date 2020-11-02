	options(scipen=999)
	project_directory <- "/lustre/scratch/jmanthey/23_camp1"
	directory_name <- "camp_filter"
	queue <- "omni"
	cluster <- "quanah"
	output_name <- "camp1"
	vcftools_options <- "--max-missing 0.7 --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 70"
	
	# read in the vcf list file
	vcf_list <- "vcf_list.txt"
	vcf <- scan("vcf_list.txt", what="character")
	# define output base name
	out_name <- sapply(strsplit(vcf, "\\.g\\."), "[[", 1)

	# make directories
	dir.create(directory_name)
	
	# write the two helper files
	write(vcf, file=paste(directory_name, "/helper7.txt", sep=""), ncolumns=1)
	write(out_name, file=paste(directory_name, "/helper8.txt", sep=""), ncolumns=1)

	# write the array script
	a.script <- paste(directory_name, "/filter_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#$ -V", file=a.script, append=T)
	write("#$ -cwd", file=a.script, append=T)
	write("#$ -S /bin/bash", file=a.script, append=T)
	write(paste("#$ -N ", "filter_vcf", sep=""), file=a.script, append=T)
	write(paste("#$ -q ", queue, sep=""), file=a.script, append=T)
	write("#$ -pe sm 1", file=a.script, append=T)
	write(paste("#$ -P ", cluster, sep=""), file=a.script, append=T)
	write("#$ -l h_rt=12:00:00", file=a.script, append=T)
	write("#$ -l h_vmem=8G", file=a.script, append=T)
	write(paste("#$ -t 1:", length(vcf), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("input_array=$( head -n${SGE_TASK_ID} helper7.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("output_array=$( head -n${SGE_TASK_ID} helper8.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	vcftools_total <- paste("vcftools --vcf ", project_directory, "/03_vcf/${input_array} ", vcftools_options, 
		" --recode --recode-INFO-all --out ", project_directory, "/03_vcf/${output_array}", sep="")
	write(vcftools_total, file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("bcftools query -f \'%POS\\t%REF\\t%ALT[\\t%GT]\\n \' ", project_directory, "/03_vcf/${output_array}.recode.vcf > ",
			project_directory, "/03_vcf/${output_array}.simple.vcf", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("bgzip ", project_directory, "/03_vcf/${output_array}.recode.vcf", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("tabix ", project_directory, "/03_vcf/${output_array}.recode.vcf.gz", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)

