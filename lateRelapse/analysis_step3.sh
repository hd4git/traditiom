### Filter for PASS including germline ###
cd ~/lateRelapse_targetedPanel/mutect2VCFs_withPON_chr

while read f; 
do 
grep "#" $f\_tumor-only_mutect2_filtered_chr.vcf | grep -vE "##contig=<ID=chr[A-Za-z0-9]*_|##contig=<ID=HLA|##contig=<ID=chrEBV" > $f\_PASSgermline.vcf
grep -v "#" $f\_tumor-only_mutect2_filtered_chr.vcf | grep -E "PASS|germline" >> $f\_PASSgermline.vcf
bgzip $f\_PASSgermline.vcf 
tabix -p vcf $f\_PASSgermline.vcf.gz
done < <(ls | awk "{split($0,a,"_"); print a[1]"_"a[2]}")

### Restrict variant calls to targeted panel ###
while read f; 
do 
bcftools view $f\.vcf.gz --regions-file ../src/oncomineCoveredByProbesHg38.bed > $f\_targeted.vcf 
bgzip $f\_targeted.vcf 
tabix -p vcf $f\_targeted.vcf.gz 
done < <(ls *vcf.gz | awk "{split($0, a, "."); print a[1]}")
