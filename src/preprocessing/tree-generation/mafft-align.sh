#!/bin/bash

clear

date 


# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=test
#SBATCH --job-name=mafft
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=scripts/preprocessing/tree-generation/jobLog_%J.out
#SBATCH --error=scripts/preprocessing/tree-generation/jobLog_%J.err
#############################################


infile=data/abund-tax-raw/asv_blanes_16s-short.fasta


module load mafft/7.402

mafft --thread 24 --reorder --auto  ${infile} > data/abund-tax-raw/blanes.pir
