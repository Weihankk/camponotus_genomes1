library(Biostrings)

library(stringr)

options(scipen=999)



# read in genome

genome <- readDNAStringSet("camponotus_reordered.FINAL.fasta")

# get contig names and lengths

genome_names <- genome@ranges@NAMES

genome_lengths <- genome@ranges@width



# sort by genome lengths

genome <- genome[order(genome_lengths, decreasing=T)]

# get names and lengths again to check

genome_names <- genome@ranges@NAMES

genome_lengths <- genome@ranges@width

head(genome_lengths)



# remake contig names to be simple

genome@ranges@NAMES <- paste("scaffold", str_pad(seq(from=1, to=length(genome_names), by=1), 4, pad="0"), sep="")

# remove contigs less than 20kb 
genome2 <- genome[genome@ranges@width >= 20000]

# write output

writeXStringSet(genome2, file="camponotus_sp_genome.fasta")
