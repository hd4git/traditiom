#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=1:mem=8gb

src=$RDS/ephemeral/barcodes/raw/X204SC20063293-Z01-F001/raw_data/
results=$RDS/ephemeral/barcodes/results
ref=$RDS/ephemeral/barcodes/ref/cloneTracker_lib_191014.fa
scripts=/rds/general/user/hdhiman/home/scripts/barcode_analysis


module load anaconda3/personal
source activate genomics
i=1

echo "QC raw data..."
mkdir $results/fastqc/raw/S$i
fastqc $src/S$i/S$i*.fq.gz  -o $results/fastqc/raw/S$i -d $results/fastqc/tmp

echo "Trim raw data..."
mkdir $results/fastqc/trimmed/S$i
mkdir $results/trimmed/S$i/
trim_galore --paired $src/S$i/S$i*_1.fq.gz $src/S$i/S$i*_2.fq.gz --phred33 -q 30 --gzip -o $results/trimmed/S$i/
fastqc $results/trimmed/S$i/S$i*val*.fq  -o $results/fastqc/trimmed/S$i -d $results/fastqc/tmp

echo "Merging fastq..."
cat $results/trimmed/S$i/S$i*.fq > $results/merged_fq/S$i\_merged.fq 
