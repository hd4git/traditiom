~/T47D_targetedPanel/T47D_platypus_joined_mutation_calls_PASS_sorted.vcf.gz


~/wgs/timon/T47D_patypus/T47D_patypus_filtered.vcf
awk '{print $1"\t"$2"\t"$2+1}' ~/wgs/timon/T47D_patypus/T47D_patypus_filtered.vcf > ~/wgs/timon/T47D_patypus/T47D_patypus_filtered.bed

vcftools --gzvcf ~/T47D_targetedPanel/T47D_platypus_joined_mutation_calls_PASS_sorted.vcf.gz --bed <(awk '{print $1"\t"$2"\t"$2+1}' ~/wgs/timon/T47D_patypus/T47D_patypus_filtered.vcf) --out ~/wgs/timon/T47D_patypus/T47D_patypus_filtered_overllap.vcf --recode --keep-INFO-all




grep -v "#" T47D_platypus_joined_mutation_calls_PASS_sorted.vcf | awk '{print $1":"$2"_"$4"/"$5"\t"$0}' > T47D_platypus_joined_mutation_calls_PASS_sorted2.vcf
awk 'NR==FNR{a[$1]=$0; next} ($1 in a) {print $0}' T47D_patypus_filtered.txt T47D_platypus_joined_mutation_calls_PASS_sorted2.vcf > T47D_platypus_filtered_ori.txt





grep -v "#" T47D_platypus_joined_mutation_calls_PASS_sorted.vcf | awk '{print $1":"$2"\t"$0}' > T47D_platypus_joined_mutation_calls_PASS_sortedAll.vcf
grep "#" ~/wgs/timon/T47D_patypus/T47D_platypus_joined_mutation_calls_PASS_sorted.vcf > T47D_platypus_filtered_oriAll.txt
awk 'NR==FNR{a[$1]=$0; next} ($1 in a) {print $0}' <(awk '{split($0, a, "_"); print a[1]}' T47D_patypus_filtered.txt) T47D_platypus_joined_mutation_calls_PASS_sortedAll.vcf | cut -f 2- >> T47D_platypus_filtered_oriAll.txt


