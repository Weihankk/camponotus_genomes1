library(Biostrings)

library(stringr)

options(scipen=999)



# read in genome

genome <- readDNAStringSet("camp.contigs.fasta")

# get contig names and lengths

genome_names <- genome@ranges@NAMES

genome_lengths <- as.numeric(sapply(strsplit(sapply(strsplit(genome_names, " "), "[[", 2), "="), "[[", 2))



# sort by genome lengths

genome <- genome[order(genome_lengths, decreasing=T)]

# get names and lengths again

genome_names <- genome@ranges@NAMES

genome_lengths <- as.numeric(sapply(strsplit(sapply(strsplit(genome_names, " "), "[[", 2), "="), "[[", 2))

head(genome_lengths)



# remake contig names to be simple

genome@ranges@NAMES <- paste("contig", str_pad(seq(from=1, to=length(genome_names), by=1), 4, pad="0"), sep="")

# remove contigs less than 10kb 
genome2 <- genome[genome@ranges@width >= 10000]

# write output

writeXStringSet(genome2, file="camponotus_reordered.fasta")



# contig sizes for juicer:

sizes <- cbind(genome@ranges@NAMES, genome_lengths)

write.table(sizes, file="camponotus_reordered.contig.sizes", sep="\t", quote=F, row.names=F, col.names=F)
