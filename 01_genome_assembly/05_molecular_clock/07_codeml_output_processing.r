
output_name <- "codeml_results_pvals_uncorrected.txt"
output_name2 <- "codeml_results_pvals_corrected.txt"
write(c("gene_number", "camponotus_un_p", "formica_un_p", "nylanderia_un_p", "lasius_un_p"), file=output_name, ncolumns=5, sep="\t")

x <- list.files(pattern="*fasta")
x2 <- list.files(pattern="*txt")
x_numbers <- as.numeric(sapply(strsplit(x, "_"), "[[", 1))

library(stats)
library(Biostrings)
library(seqinr)

for(a in 1:max(x_numbers)) {
	x_match <- match(paste(a, "_aligned_trimmed.fasta", sep=""), x) # see if fasta file exists
	if(!is.na(x_match)) { # if yes read it
		a_rep <- readDNAStringSet(paste(a, "_aligned_trimmed.fasta", sep=""))
		if(a_rep@ranges@width[1] >= 150) { # only use those that are at least 50 AAs (150 nucleotides)
			x_match <- match(paste(a, "_total_output.txt", sep=""), x2)
			if(!is.na(x_match)) { # read in codeml output
				a_results <- read.table(paste(a, "_total_output.txt", sep=""), fill = T, stringsAsFactors=F)
				a_null_lnl <- a_results[1,5] # get the null model lnL
				a_alt_lnl <- a_results[c(3,5,7,9), 5] # get the alt models' lnL
				a_LRT <- 2 * (a_alt_lnl - a_null_lnl) # calculate the LRT
				a_uncorrected_p <- pchisq(a_LRT, df=1, lower.tail=FALSE) # get p-values for the LRT values (chi-square two tail)
				a_output <- c(a, a_uncorrected_p)
				write(a_output, file=output_name, ncolumns=5, append=T, sep="\t")
			}
		}
	} 
}

# read in previous output with p-values
output <- read.table(output_name, sep="\t", stringsAsFactors=F, header=T)
# calculate number of tests = number of genes * 4 tests
number_comparisons <- nrow(output) * 4
# multiple testing correction of the p-values using Benjamini & Hochberg (1995) (fdr)
output[,2] <- p.adjust(output[,2], method="fdr", n=number_comparisons)
output[,3] <- p.adjust(output[,3], method="fdr", n=number_comparisons)
output[,4] <- p.adjust(output[,4], method="fdr", n=number_comparisons)
output[,5] <- p.adjust(output[,5], method="fdr", n=number_comparisons)


# find minimum p-value for each gene and append that column to the output
min_p <- apply(output[,2:5], 1, min)
output <- cbind(output, min_p)
plot(min_p, pch=19, cex=0.1)

write.table(output, file=output_name2, sep="\t", quote=F, col.names=T, row.names=F)



# make the 4 fold degenerate sites output directory
dir.create("_4d_output")

# remove all significant tests and any rows missing info
filtered_output <- na.omit(output)
filtered_output <- filtered_output[filtered_output$min_p > 0.05, ]

# function to identify 4 fold degenerate sites
codon_table <- read.table("codon_table.txt", header=T, stringsAsFactors=F)
determine_4d <- function(xxxx) {
	# remove codons with N
	if(length(grep("N", xxxx)) == 0) {
		# check if matches have 4D sites and return them
		if(codon_table[match(xxxx[1], codon_table[,1]),3] == "no") {
			return("")
		} else if(codon_table[match(xxxx[1], codon_table[,1]),2] == codon_table[match(xxxx[2], codon_table[,1]),2] & codon_table[match(xxxx[1], codon_table[,1]),2] == codon_table[match(xxxx[3], codon_table[,1]),2] & codon_table[match(xxxx[1], codon_table[,1]),2] == codon_table[match(xxxx[4], codon_table[,1]),2]) {
			return_object <- substr(xxxx, 3, 3)
		} else {
			return("")
		}
	}
}

# loop to
# read in multiple sequence alignments that are not under selection so as to get the four-fold degenerate sites
# for a later phylogeny
locus_results <- list()
for(a in 1:nrow(filtered_output)) {
	a_rep <- read.fasta(paste(filtered_output[a,1], "_aligned_trimmed.fasta", sep=""))
	
	# determine how many codons and put into a list for each codon
	codons <- length(a_rep$camponotus) / 3
	
	# make codon list
	codon_list <- list()
	for(b in 1:codons) {
		codon_list[[b]] <- c(paste(toupper(a_rep$camponotus[(b*3-2):(b*3)]), collapse=""),
							paste(toupper(a_rep$lasius[(b*3-2):(b*3)]), collapse=""),
							paste(toupper(a_rep$formica[(b*3-2):(b*3)]), collapse=""),
							paste(toupper(a_rep$nylanderia[(b*3-2):(b*3)]), collapse=""))
	}
	
	# use the determine_4d function to return four fold degenerate sites for this locus
	fourd_sites <- list()
	for(b in 1:length(codon_list)) {
		fourd_sites[[b]] <- determine_4d(codon_list[[b]])
	}
	# remove null results
	fourd_sites2 <- list()
	for(b in 1:length(fourd_sites)) {
		if(b == 1) { b_count <- 1 }
		 if(length(fourd_sites[[b]]) > 1) {
			fourd_sites2[[b_count]] <- fourd_sites[[b]]
			b_count <- b_count + 1
		}
	}
	# add sites to overall results list
	locus_results[[a]] <- c(paste(sapply(fourd_sites2, "[[", 1), collapse=""),
							paste(sapply(fourd_sites2, "[[", 2), collapse=""),
							paste(sapply(fourd_sites2, "[[", 3), collapse=""),
							paste(sapply(fourd_sites2, "[[", 4), collapse=""))
	
	# print progress
	if(a %% 100 == 0) {
		print(a / nrow(filtered_output))
	}
}

# concatenate results across loci
camponotus <- paste(sapply(locus_results, "[[", 1), collapse="")
lasius <- paste(sapply(locus_results, "[[", 2), collapse="")
formica <- paste(sapply(locus_results, "[[", 3), collapse="")
nylanderia <- paste(sapply(locus_results, "[[", 4), collapse="")
# total of 806844 sites

output_name <- "_total_4d_sites.fasta"
write(">camponotus", file=output_name, ncolumns=1)
write(camponotus, file=output_name, ncolumns=1, append=T)
write(">lasius", file=output_name, ncolumns=1, append=T)
write(lasius, file=output_name, ncolumns=1, append=T)
write(">formica", file=output_name, ncolumns=1, append=T)
write(formica, file=output_name, ncolumns=1, append=T)
write(">nylanderia", file=output_name, ncolumns=1, append=T)
write(nylanderia, file=output_name, ncolumns=1, append=T)








