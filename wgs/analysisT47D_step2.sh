#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:mem=100gb

module load anaconda3/personal
source activate jvCalls

ref=/rds/general/project/traditiom/live/ref/wgs/broad_bundle_hg38_20180205
genome_snp=$ref/dbsnp_146.hg38.vcf.gz
germline_file=$ref/af-only-gnomad.hg38.vcf.gz
temp=/rds/general/user/hdhiman/ephemeral/temp
project=/rds/general/user/hdhiman/ephemeral/T47D_targetedPanel
results=$project/results
raw=$project/raw
merged=$results/merged
trimmed=$results/trimmed

bam_dir=$results/finalBAMs
fasta_file=$ref/Homo_sapiens_assembly38.fasta
chr=$ref/chr.txt

platypus callVariants \
  --refFile=$fasta_file \
  --bamFiles=$bam_dir/bam_list.txt \
  --output=$project/results/T47D_platypus_joined_mutation_calls.vcf.gz \
  --nCPU=10 \
  --filterDuplicates=0 \
  --maxReads=10000000
#   --regions $chr \


