library(Biostrings)

library(stringr)

options(scipen=999)



# read in genome

genome <- readDNAStringSet("camp_sp_genome.fasta")

# read in blobtools output
blob <- read.table("camp.camp_blobplot.blobDB.table.txt", sep="\t", stringsAsFactors=F)

# determine which scaffolds are in hymenoptera or no hits
keep <- blob[blob[,6] == "Hymenoptera" | blob[,6] == "no-hit",1]

# filter the genome
genome2 <- genome[genome@ranges@NAMES %in% keep]

# write output/
writeXStringSet(genome2, file="camp_sp_genome_filtered.fasta")
