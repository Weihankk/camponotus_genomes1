# compare the ant host tree to each of the blochmannia gene trees

library(ape)
library(phytools)
library(seqinr)
options(scipen=999)

# read in species tree
species_tree <- read.nexus("camp_summed.tre")
species_tree <- midpoint.root(species_tree)
# remove single quotes from tip labels
species_tree$tip.label <- gsub("'", "", species_tree$tip.label)
# remove furthest outgroup 
species_tree  <- drop.tip(species_tree, "SRX5650044")
# rename tips to coincide with blochmannia trees
species_tree$tip.label[species_tree$tip.label == "C-028"] <- "28_vicinus"
species_tree$tip.label[species_tree$tip.label == "C-005"] <- "5_vicinus"
species_tree$tip.label[species_tree$tip.label == "C-029"] <- "29_vicinusunk"
species_tree$tip.label[species_tree$tip.label == "C-024"] <- "24_vicinusunk"
species_tree$tip.label[species_tree$tip.label == "C-056"] <- "56_herculeanus"
species_tree$tip.label[species_tree$tip.label == "C-049"] <- "49_herculeanus"
species_tree$tip.label[species_tree$tip.label == "C-050"] <- "50_herculeanus"
species_tree$tip.label[species_tree$tip.label == "C-019"] <- "19_modoc"
species_tree$tip.label[species_tree$tip.label == "C-002"] <- "2_modoc"
species_tree$tip.label[species_tree$tip.label == "C-036"] <- "36_modoc"
species_tree$tip.label[species_tree$tip.label == "C-010"] <- "10_pennunk"
species_tree$tip.label[species_tree$tip.label == "C-016"] <- "16_pennunk"
species_tree$tip.label[species_tree$tip.label == "C-018"] <- "18_pennunk"
species_tree$tip.label[species_tree$tip.label == "C-003"] <- "3_laevigatus"
species_tree$tip.label[species_tree$tip.label == "C-046"] <- "46_laevigatus"
species_tree$tip.label[species_tree$tip.label == "SRX022802"] <- "floridanus"
species_tree$tip.label[species_tree$tip.label == "C-006"] <- "6_ocreatus"
species_tree$tip.label[species_tree$tip.label == "C-039"] <- "39_vicinus"
plot(species_tree)

# identify the trees to read in and then parse out their names
all_trees <- list.files("bloch_fasta_trees", pattern="RAxML*", full.names=T)
tree_names <- substr(all_trees, 55, nchar(all_trees))
tree_names <- sapply(strsplit(tree_names, "\\.1_"), "[[", 2)
tree_names <- sapply(strsplit(tree_names, "\\.tre"), "[[", 1)

# calculate the Kuhner and Felsenstein (1994) branch length score differences between blochmannia
# genes and the ant species tree
kf94 <- c()
for(a in 1:length(all_trees)) {
	a_rep <- read.tree(all_trees[a])
	a_rep <- midpoint.root(a_rep)
	# measure the Kuhner and Felsenstein (1994) branch length score
	kf94 <- c(kf94, dist.topo(species_tree, a_rep, method="score")[1])	
}

# define the species for measuring phylogenetic distance
species_names <- c("tanaemyrmex", "camponotus")
groups <- list()
groups[[1]] <- c("28_vicinus", "5_vicinus", "39_vicinus", "29_vicinusunk", "24_vicinusunk")
groups[[2]] <- c("3_laevigatus", "46_laevigatus", "56_herculeanus", "49_herculeanus", "50_herculeanus", "19_modoc", "2_modoc", "36_modoc", "10_pennunk", "16_pennunk", "18_pennunk")
groups[[3]] <- c(groups[[1]], groups[[2]])

# loop for each tree
output <- c()
host_distances <- c()
for(a in 1:length(all_trees)) {
	a_rep <- read.tree(all_trees[a])
	a_rep <- midpoint.root(a_rep)
	# loop for each species
	tree_distances <- c()
	for(g in 1:length(groups)) {
		group_of_interest <- groups[[g]]
		# calculate MRCA distance for each species (for the ant species tree)
		if(a == 1) {
			# calculate MRCA of entire group
			total_mrca <- getMRCA(species_tree, tip = group_of_interest)
			# calculate length of edges to reach common ancestor for each individual in group
			# and then average
			edge_lengths <- list()
			group_edge_lengths <- c()
			for(b in 1:length(group_of_interest)) {
				# set edge length as zero to add to
				edge_lengths[[b]] <- 0
				# what is the number of this individual?
				b_number <- match(group_of_interest[b], species_tree$tip.label)
				
				# loop throup edge table until reaching MRCA
				while_loop <- 0
				while(while_loop != total_mrca) {
					# add edge length
					edge_lengths[[b]] <- edge_lengths[[b]] + species_tree$edge.length[species_tree$edge[,2]==b_number]
			
					# change new number to that node and update the while_loop object to the node
					b_number <- species_tree$edge[species_tree$edge[,2]==b_number,1]
					while_loop <- b_number
				}
				# add the MRCA to the nodes_needed object
				group_edge_lengths <- c(group_edge_lengths, edge_lengths[[b]])
			}

			# calc the mean of the group edge lengths
			group_edge_lengths <- mean(group_edge_lengths)
			# add to vector
			host_distances <- c(host_distances, group_edge_lengths)
		}
		
		# calculate for the blochmannia tree
		# calculate MRCA of entire group
		total_mrca <- getMRCA(a_rep, tip = group_of_interest)
		# calculate length of edges to reach common ancestor for each individual in group
		# and then average
		edge_lengths <- list()
		group_edge_lengths <- c()
		for(b in 1:length(group_of_interest)) {
			# set edge length as zero to add to
			edge_lengths[[b]] <- 0
			# what is the number of this individual?
			b_number <- match(group_of_interest[b], a_rep$tip.label)
			
			# loop throup edge table until reaching MRCA
			while_loop <- 0
			while(while_loop != total_mrca) {
				# add edge length
				edge_lengths[[b]] <- edge_lengths[[b]] + a_rep$edge.length[a_rep$edge[,2]==b_number]
		
				# change new number to that node and update the while_loop object to the node
				b_number <- a_rep$edge[a_rep$edge[,2]==b_number,1]
				while_loop <- b_number
			}
			# add the MRCA to the nodes_needed object
			group_edge_lengths <- c(group_edge_lengths, edge_lengths[[b]])
		}

		# calc the mean of the group edge lengths
		group_edge_lengths <- mean(group_edge_lengths)
		# add to vector
		tree_distances <- c(tree_distances, group_edge_lengths)
	}
	output <- rbind(output, tree_distances)
}	

# mean blochmannia edge lengths / host phylogeny edge lengths
output[,1] <- output[,1] / host_distances[1]
output[,2] <- output[,2] / host_distances[2]
output[,3] <- output[,3] / host_distances[3]

# add all outputs into a data frame
output2 <- data.frame(gene_number = as.numeric(sapply(strsplit(tree_names, "_"), "[[", 1)), gene_name = as.character(sapply(strsplit(tree_names, "-"), "[[", 2)), kf94 = kf94, bloch_tanaemyrmex_edge_length = output[,1], bloch_camponotus_edge_length = output[,2], edge_length = output[,3])


# write output
write.table(output2, file="blochmannia_tree_characteristics.txt", sep="\t", quote=F, row.names=F, col.names=T)


# read in output
x <- read.table("blochmannia_tree_characteristics.txt", stringsAsFactors=F, header=T)
# order by gene number
x <- x[order(x[,1]),]
# get absolute rate estimates from Camponotus molecular evolution rate
summary(x$edge_length) * 1.983877e-09

# read in outgroup annotations
flor <- read.table("BX248583_floridanus_annotations.tsv", stringsAsFactors=F, header=T, fill=T, row.names=NULL)
# keep only annotations that are for "CDS"
flor <- flor[flor[,2] == "CDS",]
# sort
flor <- flor[order(as.numeric(flor[,4])),]
# redo row names
rownames(flor) <- seq(from=1, to=nrow(flor), by=1)

# check that gene IDs in object x match those in the floridanus annotations
test <- c()
for(a in 1:nrow(x)) {
	if(x[a,2] == flor[x[a,1],1]) {
		test <- c(test, TRUE)
	} else {
		print(a)
		test <- c(test, FALSE)
	}
}

# if the above test matched everything, continue

# extract locations for each gene (outgroup genome coordinates)
coords <- c()
for(a in 1:nrow(x)) {
	coords <- c(coords, (as.numeric(flor[x[a,1],5]) + as.numeric(flor[x[a,1],4])) / 2)
}
# add to table
x <- cbind(x[,1:2], coords, x[,3:6])

### read in alignments to obtain sequence identity for each gene in the ingroup (excluding indels)
x_files <- list.files("bloch_fasta_trees", pattern="*_aligned_trimmed.fasta$", full.names=T)
x_files_numbers <- as.numeric(sapply(strsplit(x_files, "_"), "[[", 5))

# percent identity for ingroup for each gene
pid <- list()
for(a in 1:nrow(x)) {
	# read in fasta alignment
	a_name <- x_files[match(x[a,1], x_files_numbers)]
	a_rep <- read.fasta(a_name)
	# remove outgroups
	a_rep <- a_rep[2:17]
	# pid
	pid_rep <- list()
	for(b in 1:length(a_rep[[1]])) {
		if(length(unique(as.character(lapply(a_rep, "[[", b)))) == 1) {
			pid_rep[[b]] <- 1
		} else {
			pid_rep[[b]] <- 0
		}
	}
	pid[[a]] <- sum(unlist(pid_rep)) / length(a_rep[[1]])
}
PID <- unlist(pid)
# add to data frame
x <- cbind(x, PID)

# percent identity for tanaemyrmex for each gene
pid <- list()
for(a in 1:nrow(x)) {
	# read in fasta alignment
	a_name <- x_files[match(x[a,1], x_files_numbers)]
	a_rep <- read.fasta(a_name)
	# keep only tanaemyrmex
	a_rep <- a_rep[names(a_rep) %in% groups[[1]]]
	# pid
	pid_rep <- list()
	for(b in 1:length(a_rep[[1]])) {
		if(length(unique(as.character(lapply(a_rep, "[[", b)))) == 1) {
			pid_rep[[b]] <- 1
		} else {
			pid_rep[[b]] <- 0
		}
	}
	pid[[a]] <- sum(unlist(pid_rep)) / length(a_rep[[1]])
}
PID_tanaemyrmex <- unlist(pid)
# add to data frame
x <- cbind(x, PID_tanaemyrmex)

# percent identity for camponotus for each gene
pid <- list()
for(a in 1:nrow(x)) {
	# read in fasta alignment
	a_name <- x_files[match(x[a,1], x_files_numbers)]
	a_rep <- read.fasta(a_name)
	# keep only camponotus
	a_rep <- a_rep[names(a_rep) %in% groups[[2]]]
	# pid
	pid_rep <- list()
	for(b in 1:length(a_rep[[1]])) {
		if(length(unique(as.character(lapply(a_rep, "[[", b)))) == 1) {
			pid_rep[[b]] <- 1
		} else {
			pid_rep[[b]] <- 0
		}
	}
	pid[[a]] <- sum(unlist(pid_rep)) / length(a_rep[[1]])
}
PID_camponotus <- unlist(pid)
# add to data frame
x <- cbind(x, PID_camponotus)

# write table output
write.table(x, file="blochmannia_tree_characteristics2.txt", sep="\t", quote=F, row.names=F, col.names=T)

# 20kbp means for each stat
window_coords <- seq(from=0, to=700000, by=20000)
for(a in 1:length(window_coords)) {
	# skip the last row
	if(a != length(window_coords)) {
		a_start <- window_coords[a]
		a_end <- window_coords[a+1]
		a_middle <- mean(c(window_coords[a+1], window_coords[a]))
		a_rep <- x[x$coords >= a_start & x$coords <= a_end,]
		a_rep <- c(a_middle, apply(a_rep[,4:ncol(a_rep)], 2, mean))
		if(a == 1) { 
			x_windows <- a_rep
		} else {
			x_windows <- rbind(x_windows, a_rep)
		}
	}
}
x_windows <- as.data.frame(x_windows)










#plot
par(mar=c(2,5,1,1))
par(mfrow=c(3,1))
plot(x$coords, x$kf94, pch=19, cex=0.3, xaxt="n", xlab="", col="gray70", ylab="KF94 Distance")
lines(x_windows$V1, x_windows$kf94, lty=2, lwd=2)
plot(x$coords, x$bloch_tanaemyrmex_edge_length, pch=19, cex=0.3, ylim=c(0,50), xaxt="n", xlab="", col="lightblue", ylab="Relative Evolution Rate")
points(x$coords, x$bloch_camponotus_edge_length, pch=19, cex=0.3, col="lightcoral")
points(x$coords, x$edge_length, pch=19, cex=0.3, col="gray70")
lines(x_windows$V1, x_windows$bloch_tanaemyrmex_edge_length, col="navyblue", lwd=2)
lines(x_windows$V1, x_windows$bloch_camponotus_edge_length, col="red3", lwd=2)
lines(x_windows$V1, x_windows$edge_length, lty=2, lwd=2)
plot(x$coords, x$PID_tanaemyrmex*100, pch=19, cex=0.3, col="lightblue", ylim=c(60, 100), ylab="Gene % Identity")
points(x$coords, x$PID_camponotus*100, pch=19, cex=0.3, col="lightcoral")
points(x$coords, x$PID*100, pch=19, cex=0.3, col="gray70")
lines(x_windows$V1, x_windows$PID_tanaemyrmex*100, col="navyblue", lwd=2)
lines(x_windows$V1, x_windows$PID_camponotus*100, col="red3", lwd=2)
lines(x_windows$V1, x_windows$PID*100, lty=2, lwd=2)
