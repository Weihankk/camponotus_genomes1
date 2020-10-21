# split up the depth / coverage file for subsequent processing to make a mask for smc++ 
end=$(wc -l < scaffold_list.txt | tr -d ' ')
for i in $(seq $end); 
do  
scaffold_list=$( head -n$i scaffold_list.txt | tail -n1 );
grep $scaffold_list camp_coverage.txt > ${scaffold_list}_coverage.txt;
echo $scaffold_list;
done
