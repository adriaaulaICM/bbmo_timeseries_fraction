#!/bin/sh

# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=test
#SBATCH --job-name=fasttree
#SBATCH --time=0-01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=scripts/preprocessing/tree-generation/jobLog_%J.out
#SBATCH --error=scripts/preprocessing/tree-generation/jobLog_%J.err

#############################################

clear 

date 


infile=data/abund-tax-raw/blanes.pir

scripts/programs/fasttree -gtr -nt -gamma ${infile} > data/abund-tax-raw/asv_16s_tree.nw
