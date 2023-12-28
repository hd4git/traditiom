#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=16:mem=100gb
#PBS -N lateRelapse_targetedPanel
#PBS -J 0-48

tid=$PBS_ARRAY_INDEX
module load anaconda3/personal
source activate genomics

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

### Get depth of coverage estimates ###
gatk \
DepthOfCoverage \
-R $genome_fasta \
-L /rds/general/user/hdhiman/home/wgs/OncomineDaliaLiftOver.bed \
-O $results/$f/$f\_docInfo.txt \
-I $results/$f/$f\.sorted.markdup.recal.bam
