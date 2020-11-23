

options(scipen=999)

x <- read.table("per_window_TE_cds_GC_bp.txt", sep="\t", stringsAsFactors=F, header=T)



# define window size form file
window_size <- x[1,3] - x[1,2] + 1
# define number of windows for line plots
num_windows <- 10


# what are the unique chromosomes and their bounding areas for plotting?
total_windows <- nrow(x)
chr <- unique(x[,1])
chr_polygons <- list()
# make the plotting polygons
for(a in 1:length(chr)) {
	a1 <- rownames(x)[x[,1] == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons[[a]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 0.9), c(a1, 0.9), c(a1, 0))
}



# in separate plots
par(mfrow=c(3,1))
par(mar=c(1,5,1,1))
plot(c(-1,-1), ylim=c(0,0.4), xlim=c(1, nrow(x)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Prop. CDS")
odd <- 0
for(a in 1:length(chr_polygons)) {
	if(odd == 1) {
		polygon(chr_polygons[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
points(x$cds / window_size, pch=19, cex=0.2, col="lightblue")
# sliding windows
x_axis_rep <- seq(from=1, to=total_windows, by=num_windows)
y_axis_rep <- c()
scaffold <- c()
for(b in 1:length(x_axis_rep)) {
	# pull out the rows for the window
	b_rows <- seq(from=(x_axis_rep[b] - as.integer(num_windows / 2) + 1), to=(x_axis_rep[b] + as.integer(num_windows / 2)))
	b_rows <- x[b_rows[b_rows > 0 & b_rows <= total_windows],]
	# keep only those on the same scaffold
	b_rows <- b_rows[b_rows[,1] == x[x_axis_rep[b], 1],]
	# get the mean of the stat
	y_axis_rep <- c(y_axis_rep, mean(b_rows$cds) / window_size)
	# record the scaffold
	scaffold <- c(scaffold, x[x_axis_rep[b],1])
}
# plot each scaffold at a time (so lines don't connect between scaffolds)
for(a in 1:length(unique(scaffold))) {
	a_rep <- cbind(x_axis_rep, y_axis_rep, scaffold)
	a_rep <- a_rep[a_rep[,3] == unique(scaffold)[a], ]
	lines(a_rep[,1:2], lwd=1.2, col="navyblue")
}


plot(c(-1,-1), ylim=c(0,0.8), xlim=c(1, nrow(x)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Prop. Repetitive")
odd <- 0
for(a in 1:length(chr_polygons)) {
	if(odd == 1) {
		polygon(chr_polygons[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
points(x$TE_bp / window_size, pch=19, cex=0.2, col="lightblue")
# sliding windows
x_axis_rep <- seq(from=1, to=total_windows, by=num_windows)
y_axis_rep <- c()
scaffold <- c()
for(b in 1:length(x_axis_rep)) {
	# pull out the rows for the window
	b_rows <- seq(from=(x_axis_rep[b] - as.integer(num_windows / 2) + 1), to=(x_axis_rep[b] + as.integer(num_windows / 2)))
	b_rows <- x[b_rows[b_rows > 0 & b_rows <= total_windows],]
	# keep only those on the same scaffold
	b_rows <- b_rows[b_rows[,1] == x[x_axis_rep[b], 1],]
	# get the mean of the stat
	y_axis_rep <- c(y_axis_rep, mean(b_rows$TE_bp) / window_size)
	# record the scaffold
	scaffold <- c(scaffold, x[x_axis_rep[b],1])
}
# plot each scaffold at a time (so lines don't connect between scaffolds)
for(a in 1:length(unique(scaffold))) {
	a_rep <- cbind(x_axis_rep, y_axis_rep, scaffold)
	a_rep <- a_rep[a_rep[,3] == unique(scaffold)[a], ]
	lines(a_rep[,1:2], lwd=1.2, col="navyblue")
}


plot(c(-1,-1), ylim=c(0.25,0.5), xlim=c(1, nrow(x)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Prop. GC")
odd <- 0
for(a in 1:length(chr_polygons)) {
	if(odd == 1) {
		polygon(chr_polygons[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
points(x$gc_bp / window_size, pch=19, cex=0.2, col="lightblue")
# sliding windows
x_axis_rep <- seq(from=1, to=total_windows, by=num_windows)
y_axis_rep <- c()
scaffold <- c()
for(b in 1:length(x_axis_rep)) {
	# pull out the rows for the window
	b_rows <- seq(from=(x_axis_rep[b] - as.integer(num_windows / 2) + 1), to=(x_axis_rep[b] + as.integer(num_windows / 2)))
	b_rows <- x[b_rows[b_rows > 0 & b_rows <= total_windows],]
	# keep only those on the same scaffold
	b_rows <- b_rows[b_rows[,1] == x[x_axis_rep[b], 1],]
	# get the mean of the stat
	y_axis_rep <- c(y_axis_rep, mean(b_rows$gc_bp) / window_size)
	# record the scaffold
	scaffold <- c(scaffold, x[x_axis_rep[b],1])
}
# plot each scaffold at a time (so lines don't connect between scaffolds)
for(a in 1:length(unique(scaffold))) {
	a_rep <- cbind(x_axis_rep, y_axis_rep, scaffold)
	a_rep <- a_rep[a_rep[,3] == unique(scaffold)[a], ]
	lines(a_rep[,1:2], lwd=1.2, col="navyblue")
}






# plot stat correlations
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
plot(x$TE_bp / 100000, x$cds / 100000, pch=19, cex=0.1, ylab="CDS Proportion", xlab="Repetitive Proportion")
plot(x$gc_bp / 100000, x$TE_bp / 100000, pch=19, cex=0.1, ylab="Repetitive Proportion", xlab="GC Proportion")
plot(x$gc_bp / 100000, x$cds / 100000, pch=19, cex=0.1, ylab="CDS Proportion", xlab="GC Proportion")












