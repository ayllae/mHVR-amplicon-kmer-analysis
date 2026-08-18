#!/bin/bash
#SBATCH --job-name=kmer_cluster02
#SBATCH --partition=iob_p 
#SBATCH --constraint=Genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@email.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100gb
#SBATCH --time=7-00:00:00
#SBATCH --output=28.%j.out
#SBATCH --error=28.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL

cd /working_directory/raw_files

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate /miniconda/envs/env

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Running with $SLURM_CPUS_PER_TASK CPUs"

python 01_cluster_reads.py
python 02_cluster_visualization_and_summary.py

