cd ~/lateRelapse_targetedPanel/rds

### Add headers to filtered vcf files ###
while read f; 
do 
zgrep "#" ~/lateRelapse_targetedPanel/mutect2VCFs_withPON_chr/$f\_PASSgermline_targeted.vcf.gz  > ~/lateRelapse_targetedPanel/rds/PON_PASSgermline_chr_filtered_header/$f\.vcf
cat PON_PASSgermline_chr_filtered/$f\.vcf >> ~/lateRelapse_targetedPanel/rds/PON_PASSgermline_chr_filtered_header/$f\.vcf
done < <(ls PON_PASSgermline_chr_filtered | awk "{split($0,a,"."); print a[1]}")

### Add headers to filtered vcf files ###
cd ~/lateRelapse_targetedPanel/rds/PON_PASSgermline_chr_filtered_header
while read f; 
do 
bgzip $f
tabix -p vcf $f\.gz
done < <(ls)

### Restrict variant calls to targeted panel regions ###
while read f; 
do 
bcftools view $f\.vcf.gz --regions-file ~/wgs/oncomineCoveredByProbesHg38.bed > $f\_targeted.vcf 
bgzip $f\_targeted.vcf 
tabix -p vcf $f\_targeted.vcf.gz 
done < <(ls *.vcf.gz | awk "{split($0, a, "."); print a[1]}")

### Predict variant effects with VEP ###
source activate genomics
while read f; 
do 
vcf=$f\.vcf.gz 
ofile=$f\_vep.vcf.gz
ref=/rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205
genome_fasta=/rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205/Homo_sapiens_assembly38.fasta
vep --everything \
-i $vcf \
-o $ofile \
--vcf \
--offline \
--fasta $genome_fasta \
--compress_output bgzip \
--allele_number \
--dir_cache ~/wgs/GRCh38/vep/
done < <(ls *targeted.vcf.gz | awk "{split($0, a, "."); print a[1]}")

while read f; 
do 
zgrep -vE "germline;" $f\.vcf.gz | zgrep -vE ";germline" > $f\_filtered.vcf.gz
done < <(ls *vep.vcf.gz | awk "{split($0, a, "."); print a[1]}")

### Filter variant calls based on moderate or high effect ###
while read f; 
do 
zgrep "#" $f\.vcf.gz > $f\_filtered2.vcf
zgrep -E "MODERATE|HIGH\|" $f\.vcf.gz >> $f\_filtered2.vcf
bgzip $f\_filtered2.vcf
tabix -p vcf $f\_filtered2.vcf.gz
done < <(ls *_filtered.vcf.gz | awk "{split($0, a, "."); print a[1]}")

