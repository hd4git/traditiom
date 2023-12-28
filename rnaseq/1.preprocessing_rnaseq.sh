#!/bin/bash

src=$RDS/ephemeral/rnaseq/traditiom/raw/X204SC20120982-Z01-F001/raw_data
results=$RDS/ephemeral/rnaseq/traditiom/results
scripts=$RDS/home/scripts/rnaseq/preprocessing/traditiom

tools=$RDS/home/tools

######### Pre-processing #######

while read i; 
do 
echo "
#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=48gb

module load anaconda3/personal
source activate genomics

echo "QC raw data..."
mkdir -p $results/fastqc/raw/$i
mkdir $results/fastqc/raw/$i/tmp
fastqc $src/$i/$i\_*.fq.gz  -o $results/fastqc/raw/$i -d $results/fastqc/raw/$i/tmp

echo "Trim raw data..."
mkdir -p $results/fastqc/trimmed/$i
mkdir -p $results/fastqc/trimmed/$i/tmp
mkdir -p $results/trimmed/$i/

trim_galore --paired  $src/$i/$i\_1.fq.gz $src/$i/$i\_2.fq.gz --phred33 -q 30 --gzip -o $results/trimmed/$i/

fastqc $results/trimmed/$i/$i*val*.fq.gz  -o $results/fastqc/trimmed/$i -d $results/fastqc/trimmed/$i/tmp

"> $scripts/analysis_preprocessing_$i\.sh

chmod +x $scripts/analysis_preprocessing_$i\.sh
qsub $scripts/./analysis_preprocessing_$i\.sh

done < $RDS/home/scripts/rnaseq/samples.txt

