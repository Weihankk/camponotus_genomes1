
individuals <- c("C010_pennunk", "C016_pennunk", "C018_pennunk", "C036_modoc", "C002_modoc", "C019_modoc", "C046_laevigatus", "C003_laevigatus", "C049_herculeanus","C050_herculeanus", "C056_herculeanus", "C039_vicinus",  "C005_vicinus", "C028_vicinus", "C024_vicinusunk", "C029_vicinusunk", "C006_ocreatus")
# species
species <- unique(sapply(strsplit(individuals, "_"), "[[", 2))
species
# "pennunk"     "modoc"       "laevigatus"  "herculeanus" "vicinus"     "vicinusunk"  "ocreatus" 

# read in diversity table
div <- read.table("01_heterozygosity_per_individual.txt", header=T, stringsAsFactors=F)
# remove outgroups
div <- div[div$number <= 17, ]
# for each species get diversity
heterozygosity <- c()
for(a in 1:length(species)) {
	a_rep <- div$heterozygosity[div$species == species[a]]
	a_rep <- mean(a_rep)
	heterozygosity <- c(heterozygosity, a_rep)
}
# half the heterozygosity for the MSMC r parameter
heterozygosity <- heterozygosity / 2

directory <- "/lustre/scratch/jmanthey/23_camp1/msmc_demography"

options(scipen=999)

queue <- "omni"
cluster <- "quanah"

# determine all the input directories
directories <- c()
for(a in 1:length(individuals)) {
	directories <- c(directories, paste(directory, "/", individuals[a], "/*txt", sep=""))
	# bootstraps
	for(b in 1:10) {
		directories <- c(directories, paste(directory, "/", individuals[a], "/bootstrap_", b, "/*txt", sep=""))
	}
}

# determine the output names for each 
outputs <- substr(sapply(strsplit(directories, directory), "[[", 2), 2, nchar(sapply(strsplit(directories, directory), "[[", 2)) - 5)
outputs <- gsub("/bootstrap_", "_b", outputs)

# determine the heterozygosities to use for each 
heterozygosities <- heterozygosity[match(sapply(strsplit(outputs, "_"), "[[", 2), species)]


# write helper files
write(directories, file="helper1.txt", ncolumns=1)
write(outputs, file="helper2.txt", ncolumns=1)
write(heterozygosities, file="helper3.txt", ncolumns=1)

# write the array job
a_script <- "01_ant_msmc_array.sh"
write("#!/bin/sh", file=a_script)
write("#$ -V", file=a_script, append=T)
write("#$ -cwd", file=a_script, append=T)
write("#$ -S /bin/bash", file=a_script, append=T)
write(paste("#$ -N ", "ant_dem", sep=""), file=a_script, append=T)
write(paste("#$ -q ", queue, sep=""), file=a_script, append=T)
write("#$ -pe sm 2", file=a_script, append=T)
write(paste("#$ -P ", cluster, sep=""), file=a_script, append=T)
write("#$ -l h_rt=48:00:00", file=a_script, append=T)
write("#$ -l h_vmem=8G", file=a_script, append=T)
write(paste("#$ -t 1:", length(directories), sep=""), file=a_script, append=T)
write("", file=a_script, append=T)

write("input_array=$( head -n${SGE_TASK_ID} helper1.txt | tail -n1 )", file=a_script, append=T)
write("", file=a_script, append=T)
write("output_array=$( head -n${SGE_TASK_ID} helper2.txt | tail -n1 )", file=a_script, append=T)
write("", file=a_script, append=T)
write("heter_array=$( head -n${SGE_TASK_ID} helper3.txt | tail -n1 )", file=a_script, append=T)
write("", file=a_script, append=T)

msmc_command <- "~/msmc2_linux64bit -o ${output_array} -i 20 -t 2 -m ${heter_array} -p 1*2+20*1+1*2+1*3 ${input_array}"
write(msmc_command, file=a_script, append=T)
