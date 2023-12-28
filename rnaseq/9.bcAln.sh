#!/bin/bash

home=$RDS/ephemeral/rnaseq
ref=$RDS/ephemeral/barcodes/ref/cloneTracker_lib_191014.fa
results=$RDS/ephemeral/rnaseq/traditiom/results
scripts=$RDS/home/scripts/rnaseq/bcAln/traditiom
tools=$RDS/home/tools

######### Alignment #######
while read i; 
do 
echo "
#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=48gb

module load anaconda3/personal
source activate genomics

echo "Merging fastq..."
zcat $results/trimmed/$i/$i*.fq.gz > $results/merged_fq/$i\_merged.fq 
gzip $results/merged_fq/$i\_merged.fq 

mkdir $results/bam/
echo "Starting alignment..."
bwa mem -t 10 $ref $results/merged_fq/$i\_merged.fq.gz -o $results/bam/$i\.bc.sam
samtools flagstat $results/bam/$i\.bc.sam > $results/bam/$i\.flagstat.txt

" > $scripts/analysis_$i\.sh

chmod +x $scripts/analysis_$i\.sh
qsub $scripts/./analysis_$i\.sh

done < <(ls -ltrh $results/trimmed | awk '{print $9}' | tail -12)

#$RDS/home/scripts/rnaseq/./bcAln.sh
