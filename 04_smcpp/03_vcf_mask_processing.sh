cd /lustre/scratch/jmanthey/23_camp1/03_vcf/smcpp/

# zip the bed mask
bgzip ant_smcpp_mask.txt

# tabix the mask
tabix -p bed ant_smcpp_mask.txt.gz

# copy all vcfs to cwd
cp ../*vcf .

# concatenate vcfs that are multiple files
cat header.vcf > scaffold0001.g.vcf
grep "#" -v scaffold0001__a.g.vcf >> scaffold0001.g.vcf
grep "#" -v scaffold0001__b.g.vcf >> scaffold0001.g.vcf
grep "#" -v scaffold0001__c.g.vcf >> scaffold0001.g.vcf

cat header.vcf > scaffold0002.g.vcf
grep "#" -v scaffold0002__a.g.vcf >> scaffold0002.g.vcf
grep "#" -v scaffold0002__b.g.vcf >> scaffold0002.g.vcf

cat header.vcf > scaffold0003.g.vcf
grep "#" -v scaffold0003__a.g.vcf >> scaffold0003.g.vcf
grep "#" -v scaffold0003__b.g.vcf >> scaffold0003.g.vcf

cat header.vcf > scaffold0004.g.vcf
grep "#" -v scaffold0004__a.g.vcf >> scaffold0004.g.vcf
grep "#" -v scaffold0004__b.g.vcf >> scaffold0004.g.vcf

cat header.vcf > scaffold0005.g.vcf
grep "#" -v scaffold0005__a.g.vcf >> scaffold0005.g.vcf
grep "#" -v scaffold0005__b.g.vcf >> scaffold0005.g.vcf

cat header.vcf > scaffold0006.g.vcf
grep "#" -v scaffold0006__a.g.vcf >> scaffold0006.g.vcf
grep "#" -v scaffold0006__b.g.vcf >> scaffold0006.g.vcf

cat header.vcf > scaffold0007.g.vcf
grep "#" -v scaffold0007__a.g.vcf >> scaffold0007.g.vcf
grep "#" -v scaffold0007__b.g.vcf >> scaffold0007.g.vcf

cat header.vcf > scaffold0008.g.vcf
grep "#" -v scaffold0008__a.g.vcf >> scaffold0008.g.vcf
grep "#" -v scaffold0008__b.g.vcf >> scaffold0008.g.vcf

cat header.vcf > scaffold0009.g.vcf
grep "#" -v scaffold0009__a.g.vcf >> scaffold0009.g.vcf
grep "#" -v scaffold0009__b.g.vcf >> scaffold0009.g.vcf

cat header.vcf > scaffold0010.g.vcf
grep "#" -v scaffold0010__a.g.vcf >> scaffold0010.g.vcf
grep "#" -v scaffold0010__b.g.vcf >> scaffold0010.g.vcf

cat header.vcf > scaffold0011.g.vcf
grep "#" -v scaffold0011__a.g.vcf >> scaffold0011.g.vcf
grep "#" -v scaffold0011__b.g.vcf >> scaffold0011.g.vcf

cat header.vcf > scaffold0012.g.vcf
grep "#" -v scaffold0012__a.g.vcf >> scaffold0012.g.vcf
grep "#" -v scaffold0012__b.g.vcf >> scaffold0012.g.vcf

cat header.vcf > scaffold0013.g.vcf
grep "#" -v scaffold0013__a.g.vcf >> scaffold0013.g.vcf
grep "#" -v scaffold0013__b.g.vcf >> scaffold0013.g.vcf

