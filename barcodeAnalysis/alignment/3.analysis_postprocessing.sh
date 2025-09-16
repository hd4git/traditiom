#!/bin/bash

home=$RDS/ephemeral/barcodes
src=$RDS/ephemeral/barcode_analysis/raw/X204SC21124816-Z01-F001/raw_data
results=$RDS/ephemeral/barcode_analysis/results
ref=/rds/general/project/traditiom/live/barcode_analysis/ref/cloneTracker_lib_191014.fa
scripts=/rds/general/user/hdhiman/home/scripts/barcode_analysis5

######### Post-processing #######

while read i; 
do 
echo "
#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb

module load anaconda3/personal
source activate genomics

samtools view -bhS $results/bam/$i.bc.sam -o $results/bam/$i.bam
samtools view -bh -F 2048 -q 30 $results/bam/$i.bam | samtools sort - -o $results/bam/$i.bc.sort.uniq.bam

echo "Indexing bam..."
samtools index $results/bam/$i.bc.sort.uniq.bam

echo "Creating flagstat..."
samtools flagstat $results/bam/$i.bc.sort.uniq.bam > $results/flagstat/$i.flagstat.txt

echo "Parsing bam..."
samtools view $results/bam/$i.bc.sort.uniq.bam | awk '{split(\$12, a, \":\"); if(\$3~/^bc/ && a[3]==0){print \$0}}' | cut -f3 | LC_ALL=C  sort | LC_ALL=C uniq -c | awk '{OFS=\"\t\"}{print \$2, \$1}' | sort -k1,1nr > $results/counts_filtered/$i.counts.txt

" > $scripts/analysis_postprocessing_$i.sh

chmod +x $scripts/analysis_postprocessing_$i.sh
qsub $scripts/./analysis_postprocessing_$i.sh

done < <(tail -31 $RDS/ephemeral/barcode_analysis/raw/X204SC21124816-Z01-F001/samples.txt)
