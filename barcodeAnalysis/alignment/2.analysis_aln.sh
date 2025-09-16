#!/bin/bash

home=$RDS/ephemeral/barcodes
src=$RDS/ephemeral/barcode_analysis/raw/X204SC21124816-Z01-F001/raw_data
results=$RDS/ephemeral/barcode_analysis/results
ref=/rds/general/project/traditiom/live/barcode_analysis/ref/cloneTracker_lib_191014.fa
scripts=/rds/general/user/hdhiman/home/scripts/barcode_analysis5

######### Alignment #######
while read i; 
do 
echo "
#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb

module load anaconda3/personal
source activate genomics

echo "Merging fastq..."
#zcat $results/trimmed/$i/$i*.fq.gz > $results/merged_fq/$i\_merged.fq 
gzip $results/merged_fq/$i\_merged.fq 

mkdir $results/bam/
echo "Starting alignment..."
bwa mem -t 10 $ref $results/merged_fq/$i\_merged.fq.gz -o $results/bam/$i\.bc.sam
samtools flagstat $results/bam/$i\.bc.sam > $results/bam/$i\.flagstat.txt

" > $scripts/analysis_$i\.sh

chmod +x $scripts/analysis_$i\.sh
qsub $scripts/./analysis_$i\.sh

done < <(tail -31 $RDS/ephemeral/barcode_analysis/raw/X204SC21124816-Z01-F001/samples.txt)

