

library(Hmisc)
library(stats)


# read in RELAX output
x <- read.table("01_relax_output.txt", sep="\t", header=T, stringsAsFactors=F)
# adjust p values
p_adjust <- c()
for(a in 1:nrow(x)) {
	p_adjust <- c(p_adjust, p.adjust(x$p[a], method="BY", n=6))
}
x <- cbind(x, p_adjust)

# read in diversity measures
div <- read.table("01_het_popsizes_blochgenecount.txt", header=T, stringsAsFactors=F)

species <- c("herculeanus", "laevigatus", "modoc", "pennunk", "vicinus", "vicinusunk")

test <- c()
for(a in 1:length(species)) {
	a_rep <- x[x$species == species[a],]
	
	# get mean harmonic population size for this species
	a_div <- mean(div[div$species == species[a],7])
	
	# mean k low
	mean_k_low <- mean(a_rep$k[a_rep$p_adjust <= 0.05 & a_rep$k < 1])

	# mean k high
	mean_k_high <- mean(a_rep$k[a_rep$p_adjust <= 0.05 & a_rep$k > 1])
	
	# number significant results < 1
	num_low <- nrow(a_rep[a_rep$k < 1 & a_rep$p <= 0.05, ])
	
	# number significant results > 1
	num_high <- nrow(a_rep[a_rep$k > 1 & a_rep$p <= 0.05, ])
	
	# number significant results < 1 (pval corrected)
	num_low_corrected <- nrow(a_rep[a_rep$k < 1 & a_rep$p_adjust <= 0.05, ])

	# number significant results < 1 (pval corrected)
	num_high_corrected <- nrow(a_rep[a_rep$k > 1 & a_rep$p_adjust <= 0.05, ])
		
	test <- rbind(test, c(a_div, mean_k_low, mean_k_high, num_low, num_high, num_low_corrected, num_high_corrected))

}

test <- data.frame(harmonic_pop=test[,1], mean_k_low=test[,2], mean_k_high=test[,3], num_low=test[,4], num_high=test[,5], num_low_corrected=test[,6], num_high_corrected=test[,7])
plot(test, pch=19, cex=0.8)
rcorr(as.matrix(test))
test

par(mfrow=c(2,2))
par(mar=c(4.1,4.1,1,1))
plot(test$harmonic_pop, test$num_high, pch=19, ylab="# Blochmannia Genes w/ Intensified Selection Strength", xlab="Host Species Harmonic Mean Pop. Size (1000s)", cex.lab=0.95, cex.axis=0.8)
abline(lm(test$num_high ~ test$harmonic_pop))
summary(lm(test$num_high ~ test$harmonic_pop))
plot(test$harmonic_pop, test$num_low, pch=19, ylab="# Blochmannia Genes w/ Relaxed Selection Strength", xlab="Host Species Harmonic Mean Pop. Size (1000s)", cex.lab=0.95, cex.axis=0.8)
summary(lm(test$num_low ~ test$harmonic_pop))

plot(test$harmonic_pop, test$num_high_corrected, pch=19, ylab="# Blochmannia Genes w/ Intensified Selection Strength", xlab="Host Species Harmonic Mean Pop. Size (1000s)", cex.lab=0.95, cex.axis=0.8)
abline(lm(test$num_high_corrected ~ test$harmonic_pop))
summary(lm(test$num_high_corrected ~ test$harmonic_pop))
plot(test$harmonic_pop, test$num_low_corrected, pch=19, ylab="# Blochmannia Genes w/ Relaxed Selection Strength", xlab="Host Species Harmonic Mean Pop. Size (1000s)", cex.lab=0.95, cex.axis=0.8)
summary(lm(test$num_low_corrected ~ test$harmonic_pop))






