#combine all array output into one output file

cp 1.out > total.out

for i in {2..3708}; do
cat ${i}.out >> total.out;
done




# modify the output file with R to remove info about intervals and have just scaffold names
x <- read.table("total.out", stringsAsFactors=F)

x[,1] <- sapply(strsplit(x[,1], ":"), "[[", 1)

x[,4] <- sapply(strsplit(x[,4], ":"), "[[", 1)

write.table(x, file="total_blast.out", sep="\t", row.names=F, quote=F, col.names=F)
