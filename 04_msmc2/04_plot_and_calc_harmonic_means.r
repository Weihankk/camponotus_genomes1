individuals <- c("C010_pennunk", "C016_pennunk", "C018_pennunk", "C036_modoc", "C002_modoc", "C019_modoc", "C046_laevigatus", "C003_laevigatus", "C049_herculeanus","C050_herculeanus", "C056_herculeanus", "C039_vicinus",  "C005_vicinus", "C028_vicinus", "C024_vicinusunk", "C029_vicinusunk", "C006_ocreatus")
# species
species <- unique(sapply(strsplit(individuals, "_"), "[[", 2))
species
# "pennunk"     "modoc"       "laevigatus"  "herculeanus" "vicinus"     "vicinusunk"  "ocreatus" 


# define parameters
# gen time = 2 * age of sexual maturity in closely-related species from ant literature (3-4 years depending on estimate)
# age of sex. maturity = ~3 years
# method citation = doi: 10.1016/j.cub.2015.03.047
gen <- 6
mu <- 1.983877e-09 #* gen # needs to be per generation
# plotting numbers
min_age <- 000
max_age <- 1000000
plotting_age_range <- 1000 # years for x axis


# loop through each output directory
output <- list()
output_bootstraps <- list()
for(a in 1:length(individuals)) {
	# identify bootstraps 
	a_bootstraps <- paste("output/", individuals[a], "_b", seq(from=1, to=10, by=1), ".final.txt", sep="")
	# identify main output file
	a_main <- paste("output/", individuals[a], ".final.txt", sep="")

	
	# read in main file
	a_rep <- read.table(a_main, sep="\t", header=T)
	# rearrange main file for plotting lines
	for(d in 1:nrow(a_rep)) {
		if(d == 1) {
			a_rep2 <- rbind(c(a_rep[d,2], a_rep[d,4]), c(a_rep[d,3], a_rep[d,4]))
		} else {
			a_rep2 <- rbind(a_rep2, c(a_rep[d,2], a_rep[d,4]), c(a_rep[d,3], a_rep[d,4]))
		}
	}
	a_rep <- a_rep2
	# scale by mutation rate
	a_rep[,1] <- a_rep[,1] / mu * gen
	a_rep[,2] <- (1 / a_rep[,2]) / (2 * mu)
	# remove very young and old time frames prone to error
	a_rep <- a_rep[a_rep[,1] >= min_age & a_rep[,1] <= max_age,]
	# scale by plotting age range and pop size range
	a_rep <- a_rep / plotting_age_range
	# add to output list
	output[[a]] <- a_rep
	
	# output for each bootstrap
	output_bootstraps[[a]] <- list()
	for(b in 1:length(a_bootstraps)) {
		b_rep <- read.table(a_bootstraps[b], sep="\t", header=T)
		# rearrange main file for plotting lines
		for(d in 1:nrow(b_rep)) {
			if(d == 1) {
				b_rep2 <- rbind(c(b_rep[d,2], b_rep[d,4]), c(b_rep[d,3], b_rep[d,4]))
			} else {
				b_rep2 <- rbind(b_rep2, c(b_rep[d,2], b_rep[d,4]), c(b_rep[d,3], b_rep[d,4]))
			}
		}
		b_rep <- b_rep2
		# scale by mutation rate
		b_rep[,1] <- b_rep[,1] / mu * gen
		b_rep[,2] <- (1 / b_rep[,2]) / (2 * mu)
		# remove very young and old time frames prone to error
		b_rep <- b_rep[b_rep[,1] >= min_age & b_rep[,1] <= max_age,]
		# scale by plotting age range and pop size range
		b_rep <- b_rep / plotting_age_range
		# add to output list
		output_bootstraps[[a]][[b]] <- b_rep		
	}
}



# plot all separate
a_col <- "darkgreen"
par(mfrow=c(4,5))
par(mar=c(4.5,4.5,2,0.2))
for(a in 1:length(output)) {
	plot_name1 <- individuals[a]
	plot(c(-1,1), xlim=c(50, 900), ylim=c(0,500), pch=19, cex=0.01, log="x", xlab="kya", ylab="Pop. Size (1000s)", main="", xaxt="n", yaxt="n")
	title(main=bquote(italic(.(plot_name1))), adj=0,line=0.5, cex.main=1.5)
	axis(side=2, at=c(0, 150, 300, 450), labels=T)		
	axis(side=1, at=c(20, 50, 100, 200, 500, 900), labels=T)		
	
	
	# plot bootstraps
	for(b in 1:length(output_bootstraps[[1]])) {
		lines(output_bootstraps[[a]][[b]][,1], output_bootstraps[[a]][[b]][,2], col=a_col, lwd=0.3)
	}
	lines(output[[a]][,1], output[[a]][,2], col=a_col, lwd=3)
}






# define current pop sizes 
current_pop <- c()
for(a in 1:length(output)) {
	current_pop <- c(current_pop, output[[a]][1,2])
}

# harmonic mean of pop sizes from most recent to 500k years ago
harmonic_pop <- c()
for(a in 1:length(output)) {
	out_rep <- output[[a]]
	
	# define time series
	time_series <- seq(from=as.integer(out_rep[2,1])+1, to=500, by=1)
	# time series pops
	time_pops <- c()
	for(b in 1:length(time_series)) {
		time_pops <- c(time_pops, out_rep[time_series[b] < out_rep[,1],][1,2])
	}
	# harmonic mean of this individual
	harm_rep <- length(time_pops) / sum((1 / time_pops))
	
	# add to output element
	harmonic_pop <- c(harmonic_pop, harm_rep)
}

# read in diversity table
div <- read.table("01_heterozygosity_per_individual.txt", header=T, stringsAsFactors=F)
# remove outgroups
div <- div[div$number <= 17, ]


# define output to match with diversity table
output <- data.frame(individual=as.character(gsub("C", "C-", sapply(strsplit(individuals, "_"), "[[", 1))), current_pop=as.numeric(current_pop), harmonic_pop=as.numeric(harmonic_pop))
# rearrange
output <- output[match(div$individual, output$individual),]
# combine outputs
div <- cbind(div, output)
# remove duplicate individuals column after checking
div <- div[,-6]

#write output
write.table(div, file="01_heterozygosity_and_pop_sizes.txt", sep="\t", row.names=F, quote=F)
