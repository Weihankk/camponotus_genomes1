
# read in list of fasta sequences
# when completed, put the fasta files in the yaml directory
# run the sh script in the ncbi directory that includes the pgap.py script

x <- scan("file_list.txt", what="character")

base_directory <- "/lustre/scratch/jmanthey/23_camp1/05_annotate/"

dir.create("yaml")

for(a in 1:length(x)) {
	# write generic yaml
	generic_file <- paste("yaml/generic", a, ".yaml", sep="")
	submol_file <- paste("yaml/submol", a, ".yaml", sep="")
	submol_file2 <- paste("submol", a, ".yaml", sep="")
	y <- "report_usage: true"
	write(y, file=generic_file)
	y <- "fasta:"
	write(y, file=generic_file, append=T)
	y <- "  class: File"
	write(y, file=generic_file, append=T)
	y <- paste("  location: ", x[a], sep="")
	write(y, file=generic_file, append=T)
	y <- "submol:"
	write(y, file=generic_file, append=T)
	y <- "  class: File"
	write(y, file=generic_file, append=T)
	y <- paste("  location: ", submol_file2, sep="")
	write(y, file=generic_file, append=T)
	
	# write submol yaml
	y <- "topology: 'circular'"
	write(y, file=submol_file)
	y <- "organism:"
	write(y, file=submol_file, append=T)
	y <- "  genus_species: 'Candidatus Blochmannia'"
	write(y, file=submol_file, append=T)
	y <- "contact_info:"
	write(y, file=submol_file, append=T)
	y <- "  last_name: 'Manthey'"
	write(y, file=submol_file, append=T)
	y <- "  first_name: 'Joseph'"
	write(y, file=submol_file, append=T)
	y <- "  email: 'jdmanthey@gmail.com'"
	write(y, file=submol_file, append=T)
	y <- "  organization: 'TTU'"
	write(y, file=submol_file, append=T)
	y <- "  department: 'Biology'"
	write(y, file=submol_file, append=T)
	y <- "  phone: '999-999-9999'"
	write(y, file=submol_file, append=T)
	y <- "  fax: '999-999-9999'"
	write(y, file=submol_file, append=T)
	y <- "  street: '2901 Main Street'"
	write(y, file=submol_file, append=T)
	y <- "  city: 'Lubbock'"
	write(y, file=submol_file, append=T)
	y <- "  state: 'TX'"
	write(y, file=submol_file, append=T)
	y <- "  postal_code: '79409'"
	write(y, file=submol_file, append=T)
	y <- "  country: 'USA'"
	write(y, file=submol_file, append=T)
	y <- "authors:"
	write(y, file=submol_file, append=T)
	y <- "  -     author:"
	write(y, file=submol_file, append=T)
	y <- "            last_name: 'Manthey'"
	write(y, file=submol_file, append=T)
	y <- "            first_name: 'Joseph'"
	write(y, file=submol_file, append=T)
	y <- "            middle_initial: 'D'"
	write(y, file=submol_file, append=T)
}







