
genome_index <- read.table("camp_sp_genome_filtered.fasta.fai", sep="\t", stringsAsFactors=F)
genome_index <- genome_index[,1:2]
chromosomes <- genome_index[1:31,]
options(scipen=999)

# read repeatmasker file
x <- read.table("camponotus_sp_genome_filtered.fasta.no_dups.out", sep="", stringsAsFactors=F, fill=T, header=T)

# total proportion genome
sum(x$q_end - x$q_start + 1) / sum(genome_index[,2])

window_size <- 100000
write(c("chromosome", "start", "end", "TE_bp", "mariner_bp", "gypsy_bp", "dna-unknown_bp", "dna_bp", "ltr_bp", "nonltr_bp"), file="per_window_TE_bp.txt", sep="\t", ncolumns=10)

# loop for each chromosome
for(a in 1:nrow(chromosomes)) {
	x_rep <- x[x[,5] == chromosomes[a,1], ]
	
	a_start <- 1
	a_end <- window_size
	
	n_windows <- floor(chromosomes[a,2] / window_size)
	
	# loop for each window
	for(b in 1:n_windows) {
		b_rep <- x_rep[(x_rep[,6] >= a_start & x_rep[,6] <= a_end) | (x_rep[,7] >= a_start & x_rep[,7] <= a_end), ]
		# remove portions that go on past the window
		b_rep$q_start[b_rep$q_start < a_start] <- a_start
		b_rep$q_end[b_rep$q_end > a_end] <- a_end
		
		# which match mariner or gypsy?
		mariner_rep <- b_rep[grep("DNA/Mariner_Tc1", b_rep[,11]),]
		gypsy_rep <- b_rep[grep("LTR/Gypsy", b_rep[,11]),]
		dna_unknown_rep <- b_rep[grep("DNA/Unknown", b_rep[,11]),]
		DNA_rep <- b_rep[sapply(strsplit(b_rep[,11], "/"), "[[", 1) == "DNA",]
		LTR_rep <- b_rep[sapply(strsplit(b_rep[,11], "/"), "[[", 1) == "LTR",]
		nonLTR_rep <- b_rep[sapply(strsplit(b_rep[,11], "/"), "[[", 1) == "nonLTR",]
		
		# sum the TEs
		b_rep <- sum(b_rep$q_end - b_rep$q_start + 1)
		mariner_rep <- sum(mariner_rep$q_end - mariner_rep$q_start + 1)
		gypsy_rep <- sum(gypsy_rep$q_end - gypsy_rep$q_start + 1)
		dna_unknown_rep <- sum(dna_unknown_rep$q_end - dna_unknown_rep$q_start + 1)
		DNA_rep <- sum(DNA_rep$q_end - DNA_rep$q_start + 1)
		LTR_rep <- sum(LTR_rep$q_end - LTR_rep$q_start + 1)
		nonLTR_rep <- sum(nonLTR_rep$q_end - nonLTR_rep$q_start + 1)
		
		
		# combine and then write the output	
		b_output <- c(as.character(chromosomes[a,1]), a_start, a_end, b_rep, mariner_rep, gypsy_rep, dna_unknown_rep, DNA_rep, LTR_rep, nonLTR_rep)
		write(b_output, file="per_window_TE_bp.txt", 
			sep="\t", ncolumns=10, append=T)

		# change the counters
		a_start <- a_start + window_size
		a_end <- a_end + window_size
	}
}




