#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -N T47D_targetedPanel
#PBS -J 0-25

tid=$PBS_ARRAY_INDEX
module load anaconda3/personal
source activate genomics

ref=/rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205
genome_fasta=$ref/Homo_sapiens_assembly38.fasta
genome_snp=$ref/dbsnp_146.hg38.vcf.gz
germline_file=$ref/af-only-gnomad.hg38.vcf.gz
temp=/rds/general/user/hdhiman/ephemeral/temp
project=/rds/general/user/hdhiman/ephemeral/T47D_targetedPanel
results=$project/results
raw=$project/raw
merged=$results/merged
trimmed=$results/trimmed
library="T47D_targetedPanel"

files=$(ls $raw)
arr=($files)
f=${arr[$tid]} 

echo $f 

cd $results/$f

zcat $raw/$f/*1.fq.gz > $merged/$f\_1.fq
zcat $raw/$f/*2.fq.gz > $merged/$f\_2.fq

trim_galore --paired $merged/$f\_1.fq $merged/$f\_2.fq

fastqc $merged/$f\_1.fq $merged/$f\_2.fq 
fastqc $results/$f/$f\_1_val_1.fq $results/$f/$f\_2_val_2.fq 

bwa mem -t 10 -M $genome_fasta $results/$f/$f\_1_val_1.fq $results/$f/$f\_2_val_2.fq > $results/$f/$f\.sam

sambamba view -S -h -F "not unmapped" -f bam -t 10 $results/$f/$f\.sam > $results/$f/$f\.temp.bam
picard AddOrReplaceReadGroups I=$results/$f/$f\.temp.bam O=$results/$f/$f\.bam RGID=$f RGLB=\"$library\" RGPL=NovaSeq RGPU=unit1 RGSM=$f
sambamba sort -t 16 --tmpdir=$temp -p -m 12500000000 -o $results/$f/$f\.sorted.bam $results/$f/$f\.bam
sambamba index -t 16 $results/$f/$f\.sorted.bam
sambamba markdup -t 16 -p --tmpdir=$temp --overflow-list-size=10000000 $results/$f/$f\.sorted.bam $results/$f/$f\.sorted.markdup.bam
sambamba index -t 16 $results/$f/$f\.sorted.markdup.bam
sambamba flagstat -t 16 $results/$f/$f\.sorted.markdup.bam > $results/$f/$f\.sorted.markdup.stats.txt
gatk BaseRecalibrator -I $results/$f/$f\.sorted.markdup.bam -R $genome_fasta --known-sites $genome_snp --tmp-dir $temp -O $results/$f/$f\.sorted.markdup.recal_data.table
gatk ApplyBQSR -I $results/$f/$f\.sorted.markdup.bam -R $genome_fasta --bqsr-recal-file $results/$f/$f\.sorted.markdup.recal_data.table --tmp-dir $temp -O $results/$f/$f\.sorted.markdup.recal.bam
sambamba index -t 16 -p $results/$f/$f\.sorted.markdup.recal.bam
gatk CollectSequencingArtifactMetrics -I $results/$f/$f\.sorted.markdup.recal.bam -O $results/$f/$f\.sorted.markdup.recal.tumor_artifact --FILE_EXTENSION ".txt" -R $genome_fasta

