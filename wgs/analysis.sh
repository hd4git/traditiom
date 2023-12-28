grep "#" reWGS_MCF7_platypus_joined_mutation_calls_allsamples.vcf > reWGS_MCF7_platypus_joined_mutation_calls_allsamples_PASS.vcf
grep -v "#" reWGS_MCF7_platypus_joined_mutation_calls_allsamples.vcf | \
awk '{if($7=="PASS" && $1!~"_") print $0}' | sort -k1,1 -k2,2n >> reWGS_MCF7_platypus_joined_mutation_calls_allsamples_PASS.vcf

cat reWGS_MCF7_platypus_joined_mutation_calls_allsamples_PASS.vcf > reWGS_MCF7_platypus_joined_mutation_calls_allsamples_PASS_ori.vcf




bgzip reWGS_MCF7_platypus_joined_mutation_calls_allsamples_PASS.vcf
tabix -p vcf reWGS_MCF7_platypus_joined_mutation_calls_allsamples_PASS.vcf.gz


while read C; \
do \
bcftools view -O z \
-o chr/reWGS_MCF7_platypus_allsamples_PASS.$C\.vcf.gz \
reWGS_MCF7_platypus_joined_mutation_calls_allsamples_PASS.vcf.gz \
$C ; \
done < <(zgrep -v "#" reWGS_MCF7_platypus_joined_mutation_calls_allsamples_PASS.vcf.gz | cut -f1 | sort | uniq) 

while read f;
do
tabix -p vcf $f
done < <(ls *gz)

vep_install -a cf -s homo_sapiens -y GRCh38 -c wgs/GRCh38/vep --CONVERT

while read f;
do
vep --everything \
    -i $f\.vcf.gz \
    -o $f\_vep.vcf.gz \
    --vcf \
    --offline \
    --fasta /rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205/Homo_sapiens_assembly38.fasta \
    --compress_output bgzip \
    --allele_number \
    --dir_cache ~/wgs/GRCh38/vep/
done < <(ls *gz | awk '{split($0, a, ".vcf.gz"); print a[1]}')

##########################
grep "#" T47D_platypus_joined_mutation_calls.vcf > T47D_platypus_joined_mutation_calls_PASS.vcf
grep -v "#" T47D_platypus_joined_mutation_calls.vcf | \
awk '{if($7=="PASS" && $1!~"_") print $0}' | sort -k1,1 -k2,2n >> T47D_platypus_joined_mutation_calls_PASS.vcf

bgzip T47D_platypus_joined_mutation_calls_PASS.vcf
tabix -p vcf T47D_platypus_joined_mutation_calls_PASS.vcf.gz


while read C; \
do \
bcftools view -O z \
-o chr/T47D_platypus_allsamples_PASS.$C\.vcf.gz \
T47D_platypus_joined_mutation_calls_PASS.vcf.gz \
$C ; \
done < <(zgrep -v "#" T47D_platypus_joined_mutation_calls_PASS.vcf.gz | cut -f1 | sort | uniq) 


while read f;
do
tabix -p vcf $f
done < <(ls *gz)







