# split up the depth / coverage file for subsequent processing to make a mask for smc++ 

# need to split up file for awk, use boundaries of ten scaffolds
# find total number of lines
wc -l camp_coverage.txt 
# = 321614528

# find boundary of scaffold
grep -n -m 1 'scaffold0011' camp_coverage.txt
# = 157759262

# cut the file and make another subsetted file
head -n157759261 camp_coverage.txt > camp_coverage1.txt
tail -n163855267 camp_coverage.txt > camp_temp1.txt

# find total number of lines
wc -l camp_temp1.txt 
# = 163855267

# find boundary of scaffold
grep -n -m 1 'scaffold0021' camp_temp1.txt
# = 88738280

# cut the file and make another subsetted file
head -n88738279 camp_temp1.txt > camp_coverage2.txt
tail -n75116988 camp_temp1.txt > camp_temp2.txt

# find boundary of scaffold
grep -n -m 1 'scaffold0032' camp_temp2.txt
# = 52424565

# cut the file 
head -n52424564 camp_temp2.txt > camp_coverage3.txt

# remove temp files
rm camp_temp1.txt
rm camp_temp2.txt

# pull out each scaffold from the split up file
awk -F\t '{print>$1}' camp_coverage3.txt
awk -F\t '{print>$1}' camp_coverage2.txt
awk -F\t '{print>$1}' camp_coverage1.txt


