#!/bin/bash

home=$RDS/ephemeral/rnaseq
ref=$RDS/ephemeral/rnaseq/ref/homo_sapiens/transcriptome.idx

src=$RDS/ephemeral/rnaseq/traditiom/raw/X204SC20120982-Z01-F001/raw_data
results=$RDS/ephemeral/rnaseq/traditiom/results
scripts=$RDS/home/scripts/rnaseq/traditiom
tools=$RDS/home/tools

######### Alignment #######

while read i; 
do 
echo "
#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=124gb

module load anaconda3/personal
source activate transcriptomics

export HDF5_USE_FILE_LOCKING=FALSE

kallisto quant -t 30 -i $ref -b 100 -o $results/pseudoAlignment/$i $results/trimmed/$i/$i\_1_val_1.fq.gz $results/trimmed/$i/$i\_2_val_2.fq.gz 


"> $scripts/analysis_kallisto_$i\.sh

chmod +x $scripts/analysis_kallisto_$i\.sh
qsub $scripts/./analysis_kallisto_$i\.sh

done < $RDS/home/scripts/rnaseq/missingSamples.txt


