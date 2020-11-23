library(Biostrings)
library(seqinr)
library(stringr)

options(scipen=999)

# read in window file
window_file <- read.table("per_window_TE_cds_bp.txt", stringsAsFactors=F, header=T)

# read in genome
genome <- read.fasta("camp_sp_genome_filtered.fasta")

# genome gc and at content per window empty lists
gc <- list()
at <- list()
# loop for each window
for(a in 1:nrow(window_file)) {
	a_rep <- genome[names(genome) == window_file[a,1]][[1]][window_file[a,2]:window_file[a,3]]
	gc[[a]] <- length(a_rep[a_rep == "g" | a_rep == "c"])
	at[[a]] <- length(a_rep[a_rep == "a" | a_rep == "t"])
	if(a %% 200 == 0) {
		print(a)
	}
}
gc_bp <- unlist(gc)
at_bp <- unlist(at)

# add gc and at to data frame
window_file <- cbind(window_file, gc_bp, at_bp)

# write output
write.table(window_file, file="per_window_TE_cds_GC_bp.txt", sep="\t", quote=F, row.names=F)
