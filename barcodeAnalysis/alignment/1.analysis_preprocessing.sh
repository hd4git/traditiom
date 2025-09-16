#!/bin/bash

home=$RDS/ephemeral/barcodes
src=$RDS/ephemeral/barcode_analysis/raw/X204SC21124816-Z01-F001/raw_data
results=$RDS/ephemeral/barcode_analysis/results
ref=$RDS/ephemeral/barcodes/ref/cloneTracker_lib_191014.fa
scripts=/rds/general/user/hdhiman/home/scripts/barcode_analysis5

######### Pre-processing #######

while read i; 
do 
echo "
#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb

module load anaconda3/personal
source activate genomics

echo "QC raw data..."
mkdir -p $results/fastqc/raw/$i
mkdir -p $results/fastqc/raw/$i/tmp
#fastqc $src/$i/$i\_*.fq.gz  -o $results/fastqc/raw/$i -d $results/fastqc/raw/$i/tmp

echo "Trim raw data..."
mkdir -p $results/fastqc/trimmed/$i
mkdir -p $results/fastqc/trimmed/$i/tmp
mkdir -p $results/trimmed/$i/
gunzip $src/$i/$i\_*.fq.gz 
trim_galore --paired  $src/$i/$i\_*_1.fq $src/$i/$i\_*_2.fq --phred33 -q 30 --gzip -o $results/trimmed/$i/
fastqc -t 6 $results/trimmed/$i/$i*val*.fq.gz  -o $results/fastqc/trimmed/$i -d $results/fastqc/trimmed/$i/tmp

echo "Merging fastq..."
zcat $results/trimmed/$i/$i*.fq.gz > $results/merged_fq/$i\_merged.fq 
gzip $results/merged_fq/$i\_merged.fq 

"> $scripts/analysis_preprocessing_$i\.sh

chmod +x $scripts/analysis_preprocessing_$i\.sh
qsub $scripts/./analysis_preprocessing_$i\.sh

done < <(tail -30 $RDS/ephemeral/barcode_analysis/raw/X204SC21124816-Z01-F001/samples.txt)

