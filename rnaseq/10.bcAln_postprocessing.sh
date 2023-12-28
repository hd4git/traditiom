#!/bin/bash

home=$RDS/ephemeral/rnaseq
ref=$RDS/ephemeral/barcodes/ref/cloneTracker_lib_191014.fa
results=$RDS/ephemeral/rnaseq/traditiom/results
scripts=$RDS/home/scripts/rnaseq/bcAln/traditiom
tools=$RDS/home/tools

######### Post-processing #######

while read i; 
do 
echo "
#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=1:mem=8gb

module load anaconda3/personal
source activate genomics

samtools view -bhS $results/bam/$i.bc.sam -o $results/bam/$i.bam
samtools view -bh -F 2048 -q 30 $results/bam/$i.bam | samtools sort - -o $results/bam/$i.bc.sort.uniq.bam

echo "Indexing bam..."
samtools index $results/bam/$i.bc.sort.uniq.bam

mkdir -p $results/flagstat/
echo "Creating flagstat..."
samtools flagstat $results/bam/$i.bc.sort.uniq.bam > $results/flagstat/$i.flagstat.txt

mkdir -p $results/counts_filtered/
echo "Parsing bam..."
samtools view $results/bam/$i.bc.sort.uniq.bam | awk '{split(\$12, a, \":\"); if(\$3~/^bc/ && a[3]==0){print \$0}}' | cut -f3 | sort | uniq -c | awk '{OFS=\"\t\"}{print \$2, \$1}' | sort -k1,1nr > $results/counts_filtered/$i.counts.txt

" > $scripts/analysis_postprocessing_$i.sh

chmod +x $scripts/analysis_postprocessing_$i.sh
qsub $scripts/./analysis_postprocessing_$i.sh

done < <(ls -ltrh $results/trimmed | awk '{print $9}' | tail -12)

#$RDS/home/scripts/rnaseq/./bcAln_postprocessing.sh
