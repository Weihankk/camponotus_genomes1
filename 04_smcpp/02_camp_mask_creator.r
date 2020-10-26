options(scipen=999)

# function to keep site or not
# based on rows into apply
# keep only sites that have insufficient depth at a site to make a mask
mask <- function(xxx) {
	xxx <- xxx[xxx > 7 & xxx < 200]
	if(length(xxx) > 0) {
		return(FALSE)
	}  else {
		return(TRUE)
	}
}

x_files <- list.files(pattern="scaffold0*")
# loop for each scaffold
for(a in 1:length(x_files)) {
	print(paste("Reading in Scaffold", a))
	a_rep <- read.table(x_files[a], stringsAsFactors=F)
	print(paste("Scaffold", a, "read, processing."))
	# remove outgroups
	a_rep <- a_rep[,1:19]
	# determine whether to keep sites for mask
	keep <- apply(a_rep[,3:19], 1, mask)
	# keep only sites for mask
	a_rep <- a_rep[keep,1:2]
	
	differences <- c(2,diff(a_rep[,2]))
	row_numbers <- 1:nrow(a_rep)
	row_numbers <- row_numbers[differences > 1]
	row_numbers <- sort(c(row_numbers, row_numbers -1))
	row_numbers <- row_numbers[row_numbers != 0]
	if(length(row_numbers) %% 2 != 0) {row_numbers <- c(row_numbers, nrow(a_rep))}
	a_rep <- a_rep[row_numbers,]
	# make two objects for every other row
	first_part <- a_rep[1:(nrow(a_rep) / 2) * 2 - 1, ]
	first_part[,2] <- first_part[,2] - 1
	second_part <- a_rep[1:(nrow(a_rep) / 2) * 2, ]
	output <- cbind(first_part, second_part[,2])
	
	# write output
	if(a == 1) {
		write.table(output, file="ant_smcpp_mask.txt", sep="\t", quote=F, row.names=F, col.names=F)
	} else {
		write.table(output, file="ant_smcpp_mask.txt", sep="\t", quote=F, row.names=F, col.names=F, append=T)
	}
}




