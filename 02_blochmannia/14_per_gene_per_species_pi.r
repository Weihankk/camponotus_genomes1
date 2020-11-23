# calculate nucleotide diversity for each species for each tree

library(seqinr)
options(scipen=999)

# identify the fasta files
all_fasta <- list.files("bloch_fasta_trees", pattern="*_aligned_trimmed.fasta$", full.names=T)

# label groups
groups <- list()
groups[[1]] <- c("28_vicinus", "5_vicinus", "39_vicinus")
groups[[2]] <- c("29_vicinusunk", "24_vicinusunk")
groups[[3]] <- c("3_laevigatus", "46_laevigatus")
groups[[4]] <- c("56_herculeanus", "49_herculeanus", "50_herculeanus")
groups[[5]] <- c("19_modoc", "2_modoc", "36_modoc")
groups[[6]] <- c("10_pennunk", "16_pennunk", "18_pennunk")

# return p and q allele freqs for each site in gene
# to be used with apply across columns for each slot in an alignment
# all the sites will then be summed and divided by number of sites genotyped

p_q_columns <- function(xxx) {
	# if no variation return zero
	if(length(unique(xxx)) == 1) {
		return(c(0,0))
	} else {
		# define number of chromosomes sampled 
		n_chromosomes <- length(xxx)
		# alleles = ?
		alleles <- unique(xxx)
		# p and q freqs 
		p_q <- c(length(xxx[xxx %in% alleles[1]]), length(xxx[xxx %in% alleles[2]]))
		return(p_q)	
	}
}
	

# loop for each fasta
output <- c("gene_number", "gene_name", "n_sites", "vicinus_pi", "vicinusunk_pi", "laevigatus_pi", "herculeanus_pi", "modoc_pi", "pennunk_pi")
for(a in 1:length(all_fasta)) {
	a_rep <- read.fasta(all_fasta[a])
	# extract gene number
	a_output <- sapply(strsplit(all_fasta[a], "_"), "[[", 5)
	# extract gene name
	a_output <- c(a_output, sapply(strsplit(sapply(strsplit(all_fasta[a], "_"), "[[", 6), "gene-"), "[[", 2))
	
	# loop for each group
	for(b in 1:length(groups)) {
		# subset individuals
		b_rep <- a_rep[names(a_rep) %in% groups[[b]]]
		# row bind all individuals together
		if(length(b_rep) == 2) {
			d_rep <- rbind(b_rep[[1]], b_rep[[2]])
		} else if(length(b_rep) == 3){ 
			d_rep <- rbind(b_rep[[1]], b_rep[[2]], b_rep[[3]])
		}
		# calculate pi if there are two individuals in the species or calculate pi for each pairwise comparison of two 
		# individuals and take the mean measurement for populations with 3 individuals
		# so that all species are directly comparable for sample sizes
		if(nrow(d_rep) == 2) {
			# find the p and q frequencies for each base in the gene
			p_q_rep <- t(apply(d_rep, 2, p_q_columns))
			# based on nucleotide diversity formula from: 
			# Carlson et al. 2005: 10.1101/gr.4326505 (equation 2)
			p_freq <- p_q_rep[,1]
			q_freq <- p_q_rep[,2]
			n_chromosomes <- nrow(d_rep)
			n_sites <- ncol(d_rep)
			pi <- sum(2 * p_freq * q_freq) * (n_chromosomes / (n_chromosomes - 1)) / n_sites
		} else {
			combinations <- combn(c(1,2,3), 2)
			pi_reps <- c()
			# calculate pi for each pair of comparisons
			for(d in 1:ncol(combinations)) {
				dd_rep <- d_rep[c(combinations[1,d],combinations[2,d]),]
				# find the p and q frequencies for each base in the gene
				p_q_rep <- t(apply(dd_rep, 2, p_q_columns))
				# based on nucleotide diversity formula from: 
				# Carlson et al. 2005: 10.1101/gr.4326505 (equation 2)
				p_freq <- p_q_rep[,1]
				q_freq <- p_q_rep[,2]
				n_chromosomes <- nrow(dd_rep)
				n_sites <- ncol(dd_rep)
				pi_reps <- c(pi_reps, sum(2 * p_freq * q_freq) * (n_chromosomes / (n_chromosomes - 1)) / n_sites)
			}
			pi <- mean(pi_reps)
		}
		# find the p and q frequencies for each base in the gene
		p_q_rep <- t(apply(d_rep, 2, p_q_columns))
		# based on nucleotide diversity formula from: 
		# Carlson et al. 2005: 10.1101/gr.4326505 (equation 2)
		p_freq <- p_q_rep[,1]
		q_freq <- p_q_rep[,2]
		n_chromosomes <- nrow(d_rep)
		n_sites <- ncol(d_rep)
		pi <- sum(2 * p_freq * q_freq) * (n_chromosomes / (n_chromosomes - 1)) / n_sites
		# add to output for this gene
		if(b == 1) { a_output <- c(a_output, n_sites)}
		a_output <- c(a_output, pi)
	}

	# add gene's output to the total output
	output <- rbind(output, a_output)
}
# write output preliminary table
write.table(output, file="blochmannia_pi.txt", sep="\t", row.names=F, col.names=F, quote=F)

# read in output
x <- read.table("blochmannia_pi.txt", sep="\t", header=T, stringsAsFactors=F)
# order by gene number
x <- x[order(x[,1]),]

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
x <- cbind(x[,1:2], coords, x[,3:9])

# write table output
write.table(x, file="blochmannia_pi2.txt", sep="\t", quote=F, row.names=F, col.names=T)
# read in output
x <- read.table("blochmannia_pi2.txt", sep="\t", header=T, stringsAsFactors=F)

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
par(mfrow=c(6,1))
plot(x$coords, x$vicinus_pi, pch=19, cex=0.3, xaxt="n", xlab="", col="gray70", ylab="C. vicinus pi", ylim=c(0,0.15))
lines(x_windows$V1, x_windows$vicinus_pi, lty=2, lwd=2)
plot(x$coords, x$vicinusunk_pi, pch=19, cex=0.3, xaxt="n", xlab="", col="gray70", ylab="C. vicinus unk pi", ylim=c(0,0.15))
lines(x_windows$V1, x_windows$vicinusunk_pi, lty=2, lwd=2)
plot(x$coords, x$laevigatus_pi, pch=19, cex=0.3, xaxt="n", xlab="", col="gray70", ylab="C. laevigatus pi", ylim=c(0,0.15))
lines(x_windows$V1, x_windows$laevigatus_pi, lty=2, lwd=2)
plot(x$coords, x$herculeanus_pi, pch=19, cex=0.3, xaxt="n", xlab="", col="gray70", ylab="C. herculeanus pi", ylim=c(0,0.15))
lines(x_windows$V1, x_windows$herculeanus_pi, lty=2, lwd=2)
plot(x$coords, x$modoc_pi, pch=19, cex=0.3, xaxt="n", xlab="", col="gray70", ylab="C. modoc pi", ylim=c(0,0.15))
lines(x_windows$V1, x_windows$modoc_pi, lty=2, lwd=2)
plot(x$coords, x$pennunk_pi, pch=19, cex=0.3, xaxt="n", xlab="", col="gray70", ylab="C. penn unk pi", ylim=c(0,0.15))
lines(x_windows$V1, x_windows$pennunk_pi, lty=2, lwd=2)

