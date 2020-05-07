library(ape)
library(phytools)

# total number of individuals
total_inds <- 18

# tree files
xx <- list.files(pattern="*tre$")

# fasta files aligned
fastas <- list.files(pattern="*fasta$")

# write helper file with base name of each item
base_names <- sapply(strsplit(sapply(strsplit(xx, "bipartitions."), "[[", 2), ".tre"), "[[", 1)
write(base_names, file="helper1.txt", ncolumns=1)

# check that each tree file has a matching fasta
for(a in 1:length(base_names)) {
	a_rep <- paste(base_names[a], "_aligned_trimmed.fasta", sep="")
	if(a_rep %in% fastas == FALSE) { print(a) }
}
# none missing so continue

# loop for each tree file to annotate foreground and background for RELAX analyses
# one file per tested lineage
for(a in 1:length(xx)) {
	
	# read in tree
	a_rep <- read.tree(xx[a])

	# root the tree
	a_rep <- midpoint.root(a_rep)
	a_rep <- root(a_rep, "floridanus")
	
	
	##################### laevigatus ################################################# 
	# mrca of group
	a_ingroup <- c("3_laevigatus", "46_laevigatus")
	ingroup_node_names <- "{laevigatus}"
	
	# copy the phylogeny
	a_tree <- a_rep
	
	# blank node labels object and get number of nodes
	node_labels <- list()
	all_nodes <- unique(sort(a_rep$edge[,1]))
	for(b in 1:length(all_nodes)) {
		node_labels[[b]] <- ""
	}
	
	# loop for each node
	for(b in 1:length(all_nodes)) {
		# get descendants of node
		a_desc <- getDescendants(a_tree, all_nodes[b])
		# keep only tips
		a_desc <- a_desc[a_desc <= total_inds]
		# tip labels of these descendants
		a_desc <- a_tree$tip.label[a_desc]
		# find out if these are in ingroup
		if(length(unique(a_desc %in% a_ingroup)) == 1 & unique(a_desc %in% a_ingroup)[1] == TRUE) {
			node_labels[[b]] <- ingroup_node_names
		} else {
			node_labels[[b]] <- "{background}"
		}
	}
	# loop for each tip
	for(b in 1:length(a_tree$tip.label)) {
		if(a_tree$tip.label[b] %in% a_ingroup) {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], ingroup_node_names, sep="")
		} else {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], "{background}", sep="")
		}
	}
	# add node labels to phylogeny
	a_tree$node.label <- unlist(node_labels)

	# write tree 
	write.tree(a_tree, paste(xx[a], ".newick.laevigatus", sep=""))
	
	
	
	
	
	##################### modoc ################################################# 
	# mrca of group
	a_ingroup <- c("2_modoc", "19_modoc", "36_modoc")
	ingroup_node_names <- "{modoc}"
	
	# copy the phylogeny
	a_tree <- a_rep
	
	# blank node labels object and get number of nodes
	node_labels <- list()
	all_nodes <- unique(sort(a_rep$edge[,1]))
	for(b in 1:length(all_nodes)) {
		node_labels[[b]] <- ""
	}
	
	# loop for each node
	for(b in 1:length(all_nodes)) {
		# get descendants of node
		a_desc <- getDescendants(a_tree, all_nodes[b])
		# keep only tips
		a_desc <- a_desc[a_desc <= total_inds]
		# tip labels of these descendants
		a_desc <- a_tree$tip.label[a_desc]
		# find out if these are in ingroup
		if(length(unique(a_desc %in% a_ingroup)) == 1 & unique(a_desc %in% a_ingroup)[1] == TRUE) {
			node_labels[[b]] <- ingroup_node_names
		} else {
			node_labels[[b]] <- "{background}"
		}
	}
	# loop for each tip
	for(b in 1:length(a_tree$tip.label)) {
		if(a_tree$tip.label[b] %in% a_ingroup) {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], ingroup_node_names, sep="")
		} else {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], "{background}", sep="")
		}
	}
	# add node labels to phylogeny
	a_tree$node.label <- unlist(node_labels)
	
	# write tree 
	write.tree(a_tree, paste(xx[a], ".newick.modoc", sep=""))
	
		
	##################### herculeanus ################################################# 
	# mrca of group
	a_ingroup <- c("49_herculeanus", "50_herculeanus", "56_herculeanus")
	ingroup_node_names <- "{herculeanus}"
	
	# copy the phylogeny
	a_tree <- a_rep
	
	# blank node labels object and get number of nodes
	node_labels <- list()
	all_nodes <- unique(sort(a_rep$edge[,1]))
	for(b in 1:length(all_nodes)) {
		node_labels[[b]] <- ""
	}
	
	# loop for each node
	for(b in 1:length(all_nodes)) {
		# get descendants of node
		a_desc <- getDescendants(a_tree, all_nodes[b])
		# keep only tips
		a_desc <- a_desc[a_desc <= total_inds]
		# tip labels of these descendants
		a_desc <- a_tree$tip.label[a_desc]
		# find out if these are in ingroup
		if(length(unique(a_desc %in% a_ingroup)) == 1 & unique(a_desc %in% a_ingroup)[1] == TRUE) {
			node_labels[[b]] <- ingroup_node_names
		} else {
			node_labels[[b]] <- "{background}"
		}
	}
	# loop for each tip
	for(b in 1:length(a_tree$tip.label)) {
		if(a_tree$tip.label[b] %in% a_ingroup) {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], ingroup_node_names, sep="")
		} else {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], "{background}", sep="")
		}
	}
	# add node labels to phylogeny
	a_tree$node.label <- unlist(node_labels)
	
	# write tree 
	write.tree(a_tree, paste(xx[a], ".newick.herculeanus", sep=""))
	
	##################### penn unknown ################################################# 
	# mrca of group
	a_ingroup <- c("10_pennunk", "16_pennunk", "18_pennunk")
	ingroup_node_names <- "{pennunk}"
	
	# copy the phylogeny
	a_tree <- a_rep
	
	# blank node labels object and get number of nodes
	node_labels <- list()
	all_nodes <- unique(sort(a_rep$edge[,1]))
	for(b in 1:length(all_nodes)) {
		node_labels[[b]] <- ""
	}
	
	# loop for each node
	for(b in 1:length(all_nodes)) {
		# get descendants of node
		a_desc <- getDescendants(a_tree, all_nodes[b])
		# keep only tips
		a_desc <- a_desc[a_desc <= total_inds]
		# tip labels of these descendants
		a_desc <- a_tree$tip.label[a_desc]
		# find out if these are in ingroup
		if(length(unique(a_desc %in% a_ingroup)) == 1 & unique(a_desc %in% a_ingroup)[1] == TRUE) {
			node_labels[[b]] <- ingroup_node_names
		} else {
			node_labels[[b]] <- "{background}"
		}
	}
	# loop for each tip
	for(b in 1:length(a_tree$tip.label)) {
		if(a_tree$tip.label[b] %in% a_ingroup) {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], ingroup_node_names, sep="")
		} else {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], "{background}", sep="")
		}
	}
	# add node labels to phylogeny
	a_tree$node.label <- unlist(node_labels)
	
	# write tree 
	write.tree(a_tree, paste(xx[a], ".newick.pennunk", sep=""))

	##################### vicinus unknown ################################################# 
	# mrca of group
	a_ingroup <- c("24_vicinusunk", "29_vicinusunk")
	ingroup_node_names <- "{vicinusunk}"
	
	# copy the phylogeny
	a_tree <- a_rep
	
	# blank node labels object and get number of nodes
	node_labels <- list()
	all_nodes <- unique(sort(a_rep$edge[,1]))
	for(b in 1:length(all_nodes)) {
		node_labels[[b]] <- ""
	}
	
	# loop for each node
	for(b in 1:length(all_nodes)) {
		# get descendants of node
		a_desc <- getDescendants(a_tree, all_nodes[b])
		# keep only tips
		a_desc <- a_desc[a_desc <= total_inds]
		# tip labels of these descendants
		a_desc <- a_tree$tip.label[a_desc]
		# find out if these are in ingroup
		if(length(unique(a_desc %in% a_ingroup)) == 1 & unique(a_desc %in% a_ingroup)[1] == TRUE) {
			node_labels[[b]] <- ingroup_node_names
		} else {
			node_labels[[b]] <- "{background}"
		}
	}
	# loop for each tip
	for(b in 1:length(a_tree$tip.label)) {
		if(a_tree$tip.label[b] %in% a_ingroup) {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], ingroup_node_names, sep="")
		} else {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], "{background}", sep="")
		}
	}
	# add node labels to phylogeny
	a_tree$node.label <- unlist(node_labels)
	
	# write tree 
	write.tree(a_tree, paste(xx[a], ".newick.vicinusunk", sep=""))
	
	##################### vicinus ################################################# 
	# mrca of group
	a_ingroup <- c("5_vicinus", "28_vicinus", "39_vicinus")
	ingroup_node_names <- "{vicinus}"
	
	# copy the phylogeny
	a_tree <- a_rep
	
	# blank node labels object and get number of nodes
	node_labels <- list()
	all_nodes <- unique(sort(a_rep$edge[,1]))
	for(b in 1:length(all_nodes)) {
		node_labels[[b]] <- ""
	}
	
	# loop for each node
	for(b in 1:length(all_nodes)) {
		# get descendants of node
		a_desc <- getDescendants(a_tree, all_nodes[b])
		# keep only tips
		a_desc <- a_desc[a_desc <= total_inds]
		# tip labels of these descendants
		a_desc <- a_tree$tip.label[a_desc]
		# find out if these are in ingroup
		if(length(unique(a_desc %in% a_ingroup)) == 1 & unique(a_desc %in% a_ingroup)[1] == TRUE) {
			node_labels[[b]] <- ingroup_node_names
		} else {
			node_labels[[b]] <- "{background}"
		}
	}
	# loop for each tip
	for(b in 1:length(a_tree$tip.label)) {
		if(a_tree$tip.label[b] %in% a_ingroup) {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], ingroup_node_names, sep="")
		} else {
			a_tree$tip.label[b] <- paste(a_tree$tip.label[b], "{background}", sep="")
		}
	}
	# add node labels to phylogeny
	a_tree$node.label <- unlist(node_labels)
	
	# write tree 
	write.tree(a_tree, paste(xx[a], ".newick.vicinus", sep=""))
	
}
	
