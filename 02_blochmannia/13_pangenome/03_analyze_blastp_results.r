
options(scipen=999)
library(ape)
library(phytools)
library(RColorBrewer)
library(stats)

# read in proteins file
proteins <- read.FASTA("blochmania_AA_all.fasta", type="AA")

# read in blast results
blast <- read.table("blochmania_AA_all.blast", sep="\t", stringsAsFactors=F)
colnames(blast) <- c("query", "subject", "identity", "length", "e-value")

# read in name convsersion table
name_table <- read.table("faa_gene_conversion.txt", sep="\t", stringsAsFactors=F, quote="")

# get order of individuals in phylogeny
x_files <- list.files(pattern="^output_*")
samples <- substr(x_files, 8, nchar(x_files))
# reorder samples to phylogeny ordered top to bottom
samples <- samples[c(17, 12, 9, 5, 4, 10, 13, 15, 16, 1, 2, 3, 6, 8, 7, 11, 14)]

# subset genes to >= 50% identity
blast <- blast[blast$identity >= 50,]

# duplicate blast object
blast2 <- blast

# while loop until go through and eliminate all genes
queries <- list()
matches <- list()
counter <- 1
while(nrow(blast2) > 0) {
	# subset to first query
	a_rep <- blast2[blast2$query == unique(blast2$query)[1], ]
	# get length of this query
	a_length <- length(proteins[names(proteins) == a_rep[1,1]][[1]])
	# subset blast matches to those at least 75% of query length
	a_rep <- a_rep[a_rep$length >= 0.75 * a_length,]
	
	# add queries and matches to lists
	queries[[counter]] <- a_rep[1,1]
	matches[[counter]] <- sapply(strsplit(a_rep[,2], "__"), "[[", 1)
	
	# remove query and matches from blast2 object
	blast2 <- blast2[blast2$query != a_rep[1,1] & blast2$subject != a_rep[1,1],]
	blast2 <- blast2[blast2$query %in% a_rep$subject == FALSE,]	
	blast2 <- blast2[blast2$subject %in% a_rep$subject == FALSE,]	
	
	# add one to counter 
	counter <- counter + 1
}
queries <- unlist(queries)


# make an empty gene matrix (all zeros)
# make an empty data frame with the products as row names
for(a in 1:length(samples)) {
	if(a == 1) {
		gene_matrix <- rep(0,length(queries))
	} else {
		gene_matrix <- cbind(gene_matrix, rep(0,length(queries)))
	}
}
rownames(gene_matrix) <- queries
colnames(gene_matrix) <- samples
head(gene_matrix)

# fill in matrix with presence / absence (1 / 0)
for(a in 1:length(queries)) {
	a_rep <- colnames(gene_matrix) %in% matches[[a]]
	a_rep[a_rep == T] <- 1
	a_rep[a_rep == F] <- 0
	gene_matrix[a,] <- a_rep
}
gene_matrix <- as.data.frame(gene_matrix)

# match gene products of gene_matrix and resort
gene_products <- name_table[match(rownames(gene_matrix), name_table[,1]),2]
gene_matrix <- cbind(gene_products, gene_matrix)
gene_matrix <- gene_matrix[order(gene_products),]

# write output and manually inspect (possibly merge those that were too distantly blast related to automatically match)
# also remove hypothetical proteins
write.table(gene_matrix, file="gene_pa_matrix.txt", sep="\t", quote=F, row.names=F)

x <- read.table("gene_pa_matrix_edited.txt", sep="\t", row.names=1, quote="", stringsAsFactors=F, header=T)

# get counts of total genes present
gene_counts <- as.vector(apply(x,2,sum))

# get sample names of gene counts and edit to match format of diversity table
renamed_samples <- sapply(strsplit(colnames(x), "_"), "[[", 1)
renamed_samples <- substr(renamed_samples, 2, nchar(renamed_samples))
for(a in 1:length(renamed_samples)) {
	if(nchar(renamed_samples[a]) == 1) {
		renamed_samples[a] <- paste("C-00", renamed_samples[a], sep="")
	} else if(nchar(renamed_samples[a]) == 2) {
		renamed_samples[a] <- paste("C-0", renamed_samples[a], sep="")
	} else if(nchar(renamed_samples[a]) == 3) {
		renamed_samples[a] <- paste("C-", renamed_samples[a], sep="")
	}
}

# gene counts into data frame
gene_count <- data.frame(individual=as.character(renamed_samples), gene_count=as.numeric(gene_counts))

# add new names to x object and transpose
x <- t(x)
rownames(x) <- gene_count$individual

# read in diversity table
div <- read.table("01_heterozygosity_and_pop_sizes.txt", header=T, stringsAsFactors=F)


# merge the diversity and gene count tables
# reorder
gene_count <- gene_count[match(div$individual, gene_count$individual), ]
div <- cbind(div, gene_count$gene_count)
colnames(div)[8] <- "gene_count"
# write new table
write.table(div, file="01_het_popsizes_blochgenecount.txt", sep="\t", row.names=F, quote=F)

# read in phylogeny
tree <- read.nexus("camp_summed.tre")
# remove single quotes from tip labels
tree$tip.label <- gsub("'", "", tree$tip.label)
# midpoint root
tree <- midpoint.root(tree)
# remove outgroups
tree <- drop.tip(tree, c("SRX5650044", "SRX022802"))
# rotate tree for plotting heatmap in same phylogenetic conformation as other figure
nodes_to_rotate <- c(18,19,20,21,22,25,28,29,30,31,32,33)
for(a in nodes_to_rotate) {
	tree <- rotate(tree, a)
}

#PIC
pop_size <- setNames(div[,"harmonic_pop"], div$individual)
gene_counts <- setNames(div[,"gene_count"], div$individual)
pic_pop_size <- pic(pop_size, tree)
pic_gene_counts <- pic(gene_counts, tree)
summary(lm(pic_gene_counts ~ pic_pop_size))
plot(pic_pop_size, pic_gene_counts, ylab="Blochmannia Gene Counts (PIC)", xlab="Host Harmonic Mean Pop. Size (PIC)")

# test phylogenetic signal of gene count
signal_lambda <- phylosig(tree, method="lambda", gene_counts, test=T)
signal_k <- phylosig(tree, method="K", gene_counts, test=T)

# test phylogenetic signal for each gene
ps_per_gene <- c("gene_name", "lambda", "pval")
for(a in 1:ncol(x)) {
	# extract single gene
	a_rep <- x[,a]
	
	# continue if more than one value 
	if(length(unique(a_rep)) > 1) {
		# run test (just the lambda)
		a_signal_lambda <- phylosig(tree, method="lambda", a_rep, test=T)
		a_output <- c(colnames(x)[a], a_signal_lambda$lambda, a_signal_lambda$P)
		ps_per_gene <- rbind(ps_per_gene, a_output)
	}
}
# make data frame of output
ps_per_gene <- data.frame(gene_name=as.character(ps_per_gene[2:nrow(ps_per_gene),1]), 
							lambda=as.numeric(ps_per_gene[2:nrow(ps_per_gene),2]),
							pval=as.numeric(ps_per_gene[2:nrow(ps_per_gene),3]))
							
# add corrected p values to data frame
pval_corrected <- p.adjust(ps_per_gene$pval, method="BH")
ps_per_gene <- cbind(ps_per_gene, pval_corrected)
# write table output
write.table(ps_per_gene, file="gene_phylogenetic signal.txt", sep="\t", row.names=F, col.names=T, quote=F)
# which are significant?
sig_ps_per_gene <- ps_per_gene[ps_per_gene$pval_corrected < 0.05, ]
sig_x <- x[,colnames(x) %in% sig_ps_per_gene$gene_name]


# plot heatmap of genes 
pdf(file="phylo_heatmap.pdf", height=2.5, width=6.5)
phylo.heatmap(tree,x,standardize=F, legend=F, labels=F, colors=c(brewer.pal(11,"RdYlGn")[c(7,11)]), fsize=0.7, split=c(0.3,0.7))
dev.off()

# plot heatmap of genes with significant phylogenetic signal
pdf(file="sig_phylo_heatmap.pdf", height=2.5, width=6.5)
phylo.heatmap(tree, sig_x,standardize=F, legend=F, labels=F, colors=c(brewer.pal(11,"RdYlGn")[c(7,11)]), fsize=0.7, split=c(0.3,0.7))
dev.off()


