output_directory <- "cds_fasta_files"
output_summary <- "cds_summary_and_mapping.txt"
dir.create(output_directory)

# read in fasta files
camponotus_fasta <- scan("camponotus_cds_renamed.fasta", what="character")
camponotus_fasta <- gsub(">", "", camponotus_fasta)
lasius_fasta <- scan("lasius_cds_renamed.fasta", what="character")
lasius_fasta <- gsub(">", "", lasius_fasta)
formica_fasta <- scan("formica_cds_renamed.fasta", what="character")
formica_fasta <- gsub(">", "", formica_fasta)
nylanderia_fasta <- scan("nylanderia_cds_renamed.fasta", what="character")
nylanderia_fasta <- gsub(">", "", nylanderia_fasta)

# put each fasta into a dataframe with column 1 = name and column 2 = sequence
camponotus_cds <- data.frame(name=as.character(camponotus_fasta[1:(length(camponotus_fasta) / 2) * 2 - 1]),
							sequence=as.character(camponotus_fasta[1:(length(camponotus_fasta) / 2) * 2]))
lasius_cds <- data.frame(name=as.character(lasius_fasta[1:(length(lasius_fasta) / 2) * 2 - 1]),
							sequence=as.character(lasius_fasta[1:(length(lasius_fasta) / 2) * 2]))
formica_cds <- data.frame(name=as.character(formica_fasta[1:(length(formica_fasta) / 2) * 2 - 1]),
							sequence=as.character(formica_fasta[1:(length(formica_fasta) / 2) * 2]))
nylanderia_cds <- data.frame(name=as.character(nylanderia_fasta[1:(length(nylanderia_fasta) / 2) * 2 - 1]),
							sequence=as.character(nylanderia_fasta[1:(length(nylanderia_fasta) / 2) * 2]))
							


# read in blast files
q_camponotus_s_lasius <- read.table("q_camponotus_s_lasius.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_camponotus_s_formica <- read.table("q_camponotus_s_formica.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_camponotus_s_nylanderia <- read.table("q_camponotus_s_nylanderia.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_lasius_s_camponotus <- read.table("q_lasius_s_camponotus.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_formica_s_camponotus <- read.table("q_formica_s_camponotus.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_nylanderia_s_camponotus <- read.table("q_nylanderia_s_camponotus.blast", sep="\t", stringsAsFactors=F, quote="\"")


# write header of summary file
write(c("camponotus_gene_number", "lasius_gene_name", "formica_gene_name", "nylanderia_gene_name", "number"), file=output_summary, ncolumns=5)

# list each camponotus predicted gene
camponotus_genes <- as.character(unique(camponotus_cds[,1]))

# start counter
counter <- 1

# loop for each predicted gene
for(a in 1:length(camponotus_genes)) {
	
	# matches where the subject is the other species
	s_lasius <- q_camponotus_s_lasius[grep(camponotus_genes[a], q_camponotus_s_lasius[,1]),]
	s_formica <- q_camponotus_s_formica[grep(camponotus_genes[a], q_camponotus_s_formica[,1]),]
	s_nylanderia <- q_camponotus_s_nylanderia[grep(camponotus_genes[a], q_camponotus_s_nylanderia[,1]),]
	# keep going if there is a match for all three species
	if(nrow(s_lasius) > 0 & nrow(s_formica) > 0 & nrow(s_nylanderia) > 0) {
		# check how many unique gene matches for each search
		if(length(unique(s_lasius[,2])) == 1) {
			s_lasius <- as.character(s_lasius[1,2])
		} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
			s_lasius <- s_lasius[s_lasius[,11] == min(s_lasius[,11]),]
			s_lasius <- s_lasius[s_lasius[,12] == max(s_lasius[,12]),]
			s_lasius <- as.character(s_lasius[1,2])
		}
		if(length(unique(s_formica[,2])) == 1) {
			s_formica <- as.character(s_formica[1,2])
		} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
			s_formica <- s_formica[s_formica[,11] == min(s_formica[,11]),]
			s_formica <- s_formica[s_formica[,12] == max(s_formica[,12]),]
			s_formica <- as.character(s_formica[1,2])
		}
		if(length(unique(s_nylanderia[,2])) == 1) {
			s_nylanderia <- as.character(s_nylanderia[1,2])
		} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
			s_nylanderia <- s_nylanderia[s_nylanderia[,11] == min(s_nylanderia[,11]),]
			s_nylanderia <- s_nylanderia[s_nylanderia[,12] == max(s_nylanderia[,12]),]
			s_nylanderia <- as.character(s_nylanderia[1,2])
		}
		
		# search for the matched genes in the searches where camponotus was the subject to check for reciprocal matches
		q_lasius <- q_lasius_s_camponotus[grep(s_lasius, q_lasius_s_camponotus[,1]),]
		q_formica <- q_formica_s_camponotus[grep(s_formica, q_formica_s_camponotus[,1]),]
		q_nylanderia <- q_nylanderia_s_camponotus[grep(s_nylanderia, q_nylanderia_s_camponotus[,1]),]
		# keep going if there is a match for all three species
		if(nrow(q_lasius) > 0 & nrow(q_formica) > 0 & nrow(q_nylanderia) > 0) {
			if(length(unique(q_lasius[,2])) == 1) {
				q_lasius <- as.character(q_lasius[1,2])
			} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
				q_lasius <- q_lasius[q_lasius[,11] == min(q_lasius[,11]),]
				q_lasius <- q_lasius[q_lasius[,12] == max(q_lasius[,12]),]
				q_lasius <- as.character(q_lasius[1,2])
			}
			if(length(unique(q_formica[,2])) == 1) {
				q_formica <- as.character(q_formica[1,2])
			} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
				q_formica <- q_formica[q_formica[,11] == min(q_formica[,11]),]
				q_formica <- q_formica[q_formica[,12] == max(q_formica[,12]),]
				q_formica <- as.character(q_formica[1,2])
			}
			if(length(unique(q_nylanderia[,2])) == 1) {
				q_nylanderia <- as.character(q_nylanderia[1,2])
			} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
				q_nylanderia <- q_nylanderia[q_nylanderia[,11] == min(q_nylanderia[,11]),]
				q_nylanderia <- q_nylanderia[q_nylanderia[,12] == max(q_nylanderia[,12]),]
				q_nylanderia <- as.character(q_nylanderia[1,2])
			}
			
			# keep going if all three best matches 
			if(q_formica == camponotus_genes[a] & q_nylanderia == camponotus_genes[a] & q_lasius == camponotus_genes[a]) {
				# write info to summary file
				write(c(camponotus_genes[a], s_lasius, s_formica, s_nylanderia, counter), file=output_summary, ncolumns=5, append=T)
				
				# write the fasta files for each gene that got this far
				write(">camponotus", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1)
				write(as.character(camponotus_cds[camponotus_cds[,1] == camponotus_genes[a],2]), 
					file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(">lasius", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(as.character(lasius_cds[lasius_cds[,1] == s_lasius,2]), 
					file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(">formica", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(as.character(formica_cds[formica_cds[,1] == s_formica,2]), 
					file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(">nylanderia", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(as.character(nylanderia_cds[nylanderia_cds[,1] == s_nylanderia,2]), 
					file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)	
									
				# add to counter
				counter <- counter + 1
			}
		}	
	}
}
