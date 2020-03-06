library(Biostrings)

dir.create("unaligned_fasta")
# read in files and organize
fasta_files <- list.files(pattern="*fasta")
blast_files <- list.files(pattern="*blast")

fasta_names_list <- list()
fasta_seqs_list <- list()
blast_list <- list()
for(a in 1:18) {
	if(a != 18) {
		# read in fasta
		a_rep <- scan(fasta_files[a], what="character")
		fasta_names_list[[a]] <- a_rep[grepl(">", a_rep)]
		fasta_seqs_list[[a]] <- a_rep[grepl(">", a_rep) == FALSE]
		# read in blast
		blast_list[[a]] <- read.table(blast_files[a], sep="\t", stringsAsFactors=F)
		# take only best match for each blast hit
		blast_list[[a]] <- blast_list[[a]][match(unique(blast_list[[a]][,1]), blast_list[[a]][,1]),]
		blast_list[[a]] <- blast_list[[a]][blast_list[[a]][,11] <= 0.0001,]
	} else {
		a_rep <- readDNAStringSet("BX248583_floridanus_cds.fasta")
		a_names <- a_rep@ranges@NAMES
		a_names <- sapply(strsplit(a_names, " "), "[[", 1)
		a_rep@ranges@NAMES <- a_names
		writeXStringSet(a_rep, file="floridanus_temp_cds.fasta", width=20000)
		a_rep <- scan("floridanus_temp_cds.fasta", what="character")
		reference_fasta_names <- a_rep[grepl(">", a_rep)]
		reference_fasta_seqs <- a_rep[grepl(">", a_rep) == FALSE]
	}
}

# find the floridanus genes blasted in each individual
reference_fasta_names2 <- substr(reference_fasta_names, 2, nchar(reference_fasta_names))
reference_fasta_seqs2 <- reference_fasta_seqs
for(a in 1:length(blast_list)) {
	keep <- reference_fasta_names2 %in% blast_list[[a]][,2]
	reference_fasta_names2 <- reference_fasta_names2[keep]
	reference_fasta_seqs2 <- reference_fasta_seqs2[keep]
}
length(reference_fasta_seqs2)

# extract all sequences for each gene
for(a in 1:length(reference_fasta_names2)) {
	output_name <- paste("unaligned_fasta/", reference_fasta_names2[a], ".fasta", sep="")
	output_name <- gsub("\\[", "", output_name)
	output_name <- gsub("\\]", "", output_name)
	output_name <- gsub("=", "-", output_name)
	write(">floridanus", file=output_name)
	write(reference_fasta_seqs2[a], file=output_name, append=T)
	# add seqs from each individual
	for(b in 1:length(blast_list)) {
		b_rep <- blast_list[[b]][blast_list[[b]][,2] == reference_fasta_names2[a],]
		if(nrow(b_rep) > 1) {
			b_rep <- b_rep[b_rep[,11] == min(b_rep[,11]),]
			if(nrow(b_rep) > 1) {
				b_rep == b_rep[1,]
			}
		}
		b_rep <- paste(">", b_rep[1,1], sep="")
		b_rep <- fasta_seqs_list[[b]][fasta_names_list[[b]] == b_rep]
		b_name <- paste(">", strsplit(fasta_files[b], ".fasta")[[1]][1], sep="")
		write(b_name, file=output_name, append=T)
		write(b_rep, file=output_name, append=T)
	}
}
