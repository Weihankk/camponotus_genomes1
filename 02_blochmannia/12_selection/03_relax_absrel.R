library(rjson)
library(Hmisc)
library(ape)
library(phytools)

species <- c("herculeanus", "laevigatus", "modoc", "pennunk", "vicinus", "vicinusunk")

# relax output file
out_header <- c("species", "gene", "LRT", "p", "k")
write(out_header, file="01_relax_output.txt", ncolumns=5, sep="\t")

# loop for each species
for(a in 1:length(species)) {
	a_files <- list.files(pattern=paste("relax_", species[a], "*", sep=""))
	# loop for each gene
	for(b in 1:length(a_files)) {
		b_rep <- fromJSON(file=a_files[b])
		b_output <- c(species[a],
					  sapply(strsplit(sapply(strsplit(a_files[b], "___"), "[[", 2), ".json"), "[[", 1),
					  b_rep$`test results`$LRT,
					  b_rep$`test results`$`p-value`,
					  b_rep$`test results`$`relaxation or intensification parameter`)
		write(b_output, file="01_relax_output.txt", ncolumns=5, sep="\t", append=T)
	}
}


# ABSREL parsing
x_files <- list.files(pattern="absrel___*")

# absrel output file
out_header <- c("tips_nodes", "gene", "p")
write(out_header, file="01_absrel_output.txt", ncolumns=3, sep="\t")


# loop for each file
for(a in 1:length(x_files)) {
	a_rep <- fromJSON(file=x_files[a])
	
	# find if there are any significant results
	if(a_rep$`test results`$`positive test results` > 0) {
		# read the phylogeny with node labels
		a_tree <- read.tree(text=paste(a_rep$input$trees$`0`, ";", sep=""))
		# tests to check
		tests_to_check <- names(a_rep$tested$`0`)
		# b_significant
		bsig <- c()
		# significance values
		a_pvals <- c()
		for(b in 1:length(tests_to_check)) {
			if(eval(parse(text=paste("a_rep$`branch attributes`$`0`$`", tests_to_check[b], "`$`Corrected P-value`", sep=""))) < 0.05) {
				bsig <- c(bsig, tests_to_check[b])
				a_pvals <- c(a_pvals, eval(parse(text=paste("a_rep$`branch attributes`$`0`$`", tests_to_check[b], "`$`Corrected P-value`", sep=""))))
			}
		}
		# loop for each significant result
		for(b in 1:length(bsig)) {
			b_rep <- bsig[b] 
			# check if node or tip
			if(grepl("Node", b_rep)) {
				# root tree
				b_tree <- root(a_tree, "floridanus")
				# get descendants of node
				desc <- getDescendants(b_tree, unique(b_tree$edge[,1])[b_tree$node.label == b_rep])
				# keep only tips and not other nodes
				desc <- desc[desc <= length(b_tree$tip.label)]
				# keep results only if less than all individuals (i.e., selection not at node to outgroup)
				if(length(desc) < length(b_tree$tip.label)) {
					b_out <- c(paste(b_tree$tip.label[desc], sep="__", collapse="__"), a_pvals[b], 
								sapply(strsplit(sapply(strsplit(x_files[a], "___"), "[[", 2), ".json"), "[[", 1))
					write(b_out, file="01_absrel_output.txt", ncolumns=3, sep="\t", append=T)
				}
			} else {
				# make output the tip and stats
				b_out <- c(bsig[b], a_pvals[b], sapply(strsplit(sapply(strsplit(x_files[a], "___"), "[[", 2), ".json"), "[[", 1))
				write(b_out, file="01_absrel_output.txt", ncolumns=3, sep="\t", append=T)
			}
		}
	}
}





# read in RELAX output
x <- read.table("01_relax_output.txt", sep="\t", header=T, stringsAsFactors=F)


# read in diversity measures
div <- read.table("01_heterozygosity_per_individual.txt", header=T, stringsAsFactors=F)

species <- c("herculeanus", "laevigatus", "modoc", "pennunk", "vicinus", "vicinusunk")

test <- c()
for(a in 1:length(species)) {
	a_rep <- x[x$species == species[a],]
	
	a_div <- mean(div[div$species == species[a],4])
	
	# mean k
	mean_k <- mean(a_rep$k)
	# mean significant k
	a_rep <- a_rep[a_rep$p < 0.05, ]
	mean_sig_k <- mean(a_rep$k)
	# number significant results < 1
	num_low <- nrow(a_rep[a_rep$k < 1, ])
	# number significant results > 1
	num_high <- nrow(a_rep[a_rep$k > 1, ])
	
	test <- rbind(test, c(a_div, mean_k, mean_sig_k, num_low, num_high))

}

test <- data.frame(diversity=test[,1], meank=test[,2], mean_sig_k=test[,3], num_low=test[,4], num_high=test[,5])
plot(test, pch=19, cex=0.8)
rcorr(as.matrix(test))

summary(lm(test$num_high ~ test$diversity))
par(mar=c(4.1,4.1,1,1))
plot(test$diversity, test$num_high, pch=19, ylab="# Genes w/ Intensified Selection Strength", xlab="Species Mean Observed Heterozygosity", cex.lab=0.95, cex.axis=0.8)
abline(lm(test$num_high ~ test$diversity))

