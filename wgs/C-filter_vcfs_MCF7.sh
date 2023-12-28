#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=80:00:00
#SBATCH --job-name=filter_vcfs
#SBATCH --mem-per-cpu=15000MB
#SBATCH -p smp
#SBATCH --array=0-24

set -eo pipefail

# load main conda env
module load R
module load anaconda
source activate r_env_new

TMPDIR=$(mktemp -d /rds/general/user/hdhiman/ephemeral/tmp/"$USER"__"$SLURM_JOBID"-XXXX)

#Start the job here.
Rscript ./C-filter_vcfs_MCF7.R ${SLURM_ARRAY_TASK_ID}

echo "Tmp size is (after job run): $(du -sh "$TMPDIR")"
rm -rf "$TMPDIR"