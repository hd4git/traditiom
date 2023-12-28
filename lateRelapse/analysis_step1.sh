#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=16:mem=100gb
#PBS -N lateRelapse_targetedPanel
#PBS -J 0-48

tid=$PBS_ARRAY_INDEX
module load anaconda3/personal
source activate genomics

### Set paths ###
ref=/rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205
genome_fasta=$ref/Homo_sapiens_assembly38.fasta
genome_snp=$ref/dbsnp_146.hg38.vcf.gz
germline_file=$ref/af-only-gnomad.hg38.vcf.gz
project=/rds/general/user/hdhiman/ephemeral/lateRelapse_targetedPanel
results=$project/results2
raw=$project/X204SC22120279-Z01-F001/01.RawData
trimmed=$results/trimmed
library="lateRelapse_targetedPanel"

files=$(ls $raw)
arr=($files)
f=${arr[$tid]} 

temp=$results/$f/temp
mkdir -p $results/$f
cd $results/$f

echo $f 
### Pre-processing raw reads ###
trim_galore --paired $raw/$f/*_1.fq.gz $raw/$f/*_2.fq.gz

fastqc $raw/$f/*_1.fq.gz $raw/$f/*_2.fq.gz 
fastqc $results/$f/*_1_val_1.fq.gz $results/$f/*_2_val_2.fq.gz 

### Alignment and post-processing ###
bwa mem -t 16 -M $genome_fasta $results/$f/*_1_val_1.fq.gz $results/$f/*_2_val_2.fq.gz > $results/$f/$f\.sam
samtools flagstat $results/$f/$f\.sam > $results/$f/$f\_flagstat.txt
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

### Variant calling and filtering calls ###
gatk Mutect2 -R $genome_fasta \
-I $results/$f/$f\.sorted.markdup.recal.bam \
-tumor $f \
--germline-resource "$germline_file" \
--panel-of-normals $ref/1000g_pon.hg38.vcf.gz \
--af-of-alleles-not-in-resource 0.001 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
--native-pair-hmm-threads 12 --tmp-dir $temp \
-O $results/$f/$f\_tumor-only_mutect2_unfiltered.vcf \
--bamout $results/$f/$f\_tumor-only_mutect2_unfiltered.bam \
--f1r2-tar-gz $results/$f/$f\.f1r2.tar.gz

gatk --java-options "-Xmx64G" LearnReadOrientationModel \
-I $results/$f/$f\.f1r2.tar.gz \
-O $results/$f/$f\.read-orientation-model.tar.gz \
--tmp-dir $temp 

gatk FilterMutectCalls -R $genome_fasta \
-V $results/$f/$f\_tumor-only_mutect2_unfiltered.vcf \
-O $results/$f/$f\_tumor-only_mutect2_filtered.vcf \
--ob-priors $results/$f/$f\.read-orientation-model.tar.gz \
--tmp-dir $temp 

### Restrict variant calls to chromosomes 1-22, X and Y ###
bgzip $results/mutect2VCFs_withPON/$f\_tumor-only_mutect2_filtered.vcf > $results/mutect2VCFs_withPON_chr/$f\_tumor-only_mutect2_filtered.vcf.gz
tabix -p vcf $results/mutect2VCFs_withPON/$f\_tumor-only_mutect2_filtered.vcf.gz

bcftools view $results/mutect2VCFs_withPON/$f\_tumor-only_mutect2_filtered.vcf.gz --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrM,chrX,chrY > $results/mutect2VCFs_withPON_chr/$f\_tumor-only_mutect2_filtered_chr.vcf

bgzip $results/mutect2VCFs_withPON_chr/$f\_tumor-only_mutect2_filtered_chr.vcf > $results/mutect2VCFs_withPON_chr/$f\_tumor-only_mutect2_filtered_chr.vcf.gz
tabix -p vcf  $results/mutect2VCFs_withPON_chr/$f\_tumor-only_mutect2_filtered_chr.vcf.gz

