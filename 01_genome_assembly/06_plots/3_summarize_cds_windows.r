
x <- read.table("camp_round2_cds.gff", stringsAsFactors=F)

window_file <- read.table("per_window_TE_bp.txt", stringsAsFactors=F, header=T)

# cds bp output column
cds <- list()

# loop for each window
for(a in 1:nrow(window_file)) {
	a_rep <- x[x[,1] == window_file[a,1], ]
	if(!is.null(nrow(a_rep))) {
		b_rep <- a_rep[(a_rep[,4] >= window_file[a,2] & a_rep[,4] <= window_file[a,3]) | (a_rep[,5] >= window_file[a,2] & a_rep[,5] <= window_file[a,3]), ]
		if(!is.null(nrow(b_rep))) {
			# remove portions that go on past the window
			b_rep[,4][b_rep[,4] < window_file[a,2]] <- window_file[a,2]
			b_rep[,5][b_rep[,5] > window_file[a,3]] <- window_file[a,3]
			# sum the bp in the region
			cds[[a]] <- sum(b_rep[,5] - b_rep[,4] + 1)
		} else {
			cds[[a]] <- 0
		}
	} else {
		cds[[a]] <- 0
	}
}
cds <- unlist(cds)
# add cds to data frame
window_file <- cbind(window_file, cds)
# write output
write.table(window_file, file="per_window_TE_cds_bp.txt", sep="\t", quote=F, row.names=F)



