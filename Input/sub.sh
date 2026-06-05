#!/bin/bash
#SBATCH --job-name=kmer_cluster02
#SBATCH --partition=iob_p 
#SBATCH --constraint=Genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ayllae@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100gb
#SBATCH --time=7-00:00:00
#SBATCH --output=28.%j.out
#SBATCH --error=28.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL

cd /scratch/asv10825/mHVR_June_2026/

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate /home/asv10825/miniconda/envs/kmer_env


# IMPORTANT: Avoid thread oversubscription
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Running with $SLURM_CPUS_PER_TASK CPUs"


python 02_cluster_visualization_and_summary.py

