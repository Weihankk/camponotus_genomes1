
options(scipen=999)
library(ape)

# list of pgap annotation output directories
x_files <- list.files(pattern="^output_*")
samples <- substr(x_files, 8, nchar(x_files))

# loop for each file to concatenate the fasta files and then the faa files

for(a in 1:length(x_files)) {
	# read in fasta and add it to the total fasta file
	a_fasta <- read.FASTA(list.files(x_files[a], pattern="*fasta", full.names=T))
	if(a == 1) { total_fasta <- a_fasta } else { total_fasta <- c(total_fasta, a_fasta)}
	
	# read in AA faa file
	a_faa <- read.FASTA(paste(x_files[a], "/annot.faa", sep=""), type="AA")
	
	# combine new names and annotations into table
	if(a == 1) {
		faa_gene_att <- cbind(paste(samples[a], seq(from=1, to=length(names(a_faa)), by=1), sep="__"), substr(names(a_faa), 26, nchar(names(a_faa)) - 25))
	} else {
		faa_gene_att <- rbind(faa_gene_att, cbind(paste(samples[a], seq(from=1, to=length(names(a_faa)), by=1), sep="__"), substr(names(a_faa), 26, nchar(names(a_faa)) - 25)))
	}
	
	# rename faa seqs
	names(a_faa) <- paste(samples[a], seq(from=1, to=length(names(a_faa)), by=1), sep="__")
	
	# add to final object
	if(a == 1) { total_faa <- a_faa } else { total_faa <- c(total_faa, a_faa)}
	
}

# write output
write.FASTA(total_fasta, file="blochmannia_sequences_all.fasta")
write.FASTA(total_faa, file="blochmania_AA_all.fasta")



