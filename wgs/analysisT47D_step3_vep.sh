#!/bin/bash

cd /rds/general/user/hdhiman/ephemeral/T47D_targetedPanel/data/results/final/chr
while read f;
do
echo "
#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=64gb

module load anaconda3/personal
source activate genomics

ref=/rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205
genome_fasta=$ref/Homo_sapiens_assembly38.fasta

vep --everything \
    -i /rds/general/user/hdhiman/ephemeral/T47D_targetedPanel/data/results/final/chr/$f.vcf.gz \
    -o /rds/general/user/hdhiman/ephemeral/T47D_targetedPanel/data/results/final/vep/$f\_vep.vcf.gz \
    --vcf \
    --offline \
    --fasta /rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205/Homo_sapiens_assembly38.fasta \
    --compress_output bgzip \
    --allele_number \
    --dir_cache ~/wgs/GRCh38/vep/
" > ~/T47D_targetedPanel/scripts/vep_$f\.sh

chmod +x ~/T47D_targetedPanel/scripts/vep_$f\.sh
cd ~/T47D_targetedPanel//logsVEP
qsub ~/T47D_targetedPanel/scripts/vep_$f\.sh

done < <(ls *gz | awk '{split($0, a, ".vcf.gz"); print a[1]}')

