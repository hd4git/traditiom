#!/bin/bash

while read f; 
do 
echo "
#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:mem=48gb

module load anaconda3/personal
source activate genomics

ref=/rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205
genome_fasta=$ref/Homo_sapiens_assembly38.fasta
genome_snp=$ref/dbsnp_146.hg38.vcf.gz
germline_file=$ref/af-only-gnomad.hg38.vcf.gz
project=/rds/general/user/hdhiman/ephemeral/T47D_targetedPanel
temp=/rds/general/user/hdhiman/ephemeral/temp
raw=$project/raw
merged=$project/merged
trimmed=$project/trimmed
library=\"T47D_targetedPanel\"

zcat $f/*1.fq.gz > $merged/$f\_1.fq
zcat $f/*2.fq.gz > $merged/$f\_2.fqz

gzip $merged/$f\_1.fq
gzip $merged/$f\_2.fq

trim_galore --paired $merged/$f\_1.fq.gz $merged/$f\_2.fq.gz

bwa mem -t 10 -M $genome_fasta $f\_1_val_1.fq $f\_2_val_2.fq > $f\.sam
sambamba view -S -h -F \"not unmapped\" -f bam -t 10 $f\.sam > $f\.temp.bam
picard AddOrReplaceReadGroups I=$f\.temp.bam O=$f\.bam RGID=$f RGLB=\"$library\" RGPL=NovaSeq RGPU=unit1 RGSM=$f
sambamba sort -t 16 --tmpdir=$temp -p -m 12500000000 -o $f\.sorted.bam $f\.bam
sambamba index -t 16 $f\.sorted.bam
sambamba markdup -t 16 -p --tmpdir=$temp --overflow-list-size=10000000 $f\.sorted.bam $f\.sorted.markdup.bam
sambamba index -t 16 $f\.sorted.markdup.bam
sambamba flagstat -t 16 $f\.sorted.markdup.bam > $f\.sorted.markdup.stats.txt
gatk BaseRecalibrator -I f\.sorted.markdup.bam -R $genome_fasta --known-sites $genome_snp --tmp-dir $temp -O $f\.sorted.markdup.recal_data.table
gatk ApplyBQSR -I $f\.sorted.markdup.bam -R $genome_fasta --bqsr-recal-file $f\.sorted.markdup.recal_data.table --tmp-dir temp -O $f\.sorted.markdup.recal.bam
sambamba index -t 16 -p $f\.sorted.markdup.recal.bam
gatk CollectSequencingArtifactMetrics -I $f\.sorted.markdup.recal.bam -O $f\.sorted.markdup.recal.tumor_artifact --FILE_EXTENSION \".txt\" -R $genome_fasta

" > ~/barcode_analysis/scripts/processing/gzip_$f\.sh

chmod +x ~/barcode_analysis/scripts/processing/gzip_$f\.sh
qsub ~/barcode_analysis/scripts/processing/./gzip_$f\.sh

done < <(ls -ltrh raw/ | awk '{print $9}' | sort | tail -26)

