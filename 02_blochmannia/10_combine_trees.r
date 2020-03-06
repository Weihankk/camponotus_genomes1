options(scipen=999)

# list all the files in the trees directory
x_files <- list.files(pattern="*tre$", full.names=T)

# trees into one file
tree_list <- list()
for(a in 1:length(x_files)) {
	tree_list[[a]] <- scan(x_files[a], what="character")
}
tree_list <- unlist(tree_list)
write(tree_list, file="camp_bloch_genes.trees", ncolumns=1)
