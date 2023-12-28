#!/bin/bash

home=$RDS/ephemeral/barcodes
src=$RDS/ephemeral/barcodes/raw4/X204SC21010648-Z01-F001
results=$RDS/ephemeral/barcodes/results_batch4
ref=$RDS/ephemeral/barcodes/ref/cloneTracker_lib_191014.fa
scripts=/rds/general/user/hdhiman/home/scripts/barcode_analysis_additional

######### Pre-processing #######

while read i; 
do 
echo "
#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=2:mem=64gb

module load anaconda3/personal
source activate genomics

echo "QC raw data..."
mkdir -p $results/fastqc/raw/$i
mkdir $results/fastqc/$i/tmp
fastqc $src/$i/$i\_*.fq.gz  -o $results/fastqc/raw/$i -d $results/fastqc/$i/tmp

echo "Trim raw data..."
mkdir -p $results/fastqc/trimmed/$i
mkdir -p $results/trimmed/$i/
gunzip $src/$i/$i\_*.fq.gz 
trim_galore --paired  $src/$i/$i\_*_1.fq $src/$i/$i\_*_2.fq --phred33 -q 30 --gzip -o $results/trimmed/$i/
fastqc $results/trimmed/$i/$i*val*.fq.gz  -o $results/fastqc/trimmed/$i -d $results/fastqc/$i/tmp

echo "Merging fastq..."
zcat $results/trimmed/$i/$i*.fq.gz > $results/merged_fq/$i\_merged.fq 
gzip $results/merged_fq/$i\_merged.fq 

"> $scripts/analysis_preprocessing_$i\.sh
chmod +x $scripts/analysis_preprocessing_$i\.sh
qsub $scripts/./analysis_preprocessing_$i\.sh
done < /rds/general/user/hdhiman/home/scripts/additional.txt

