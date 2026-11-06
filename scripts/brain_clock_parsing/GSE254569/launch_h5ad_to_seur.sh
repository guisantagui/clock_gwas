#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=GSE254569_2Seur
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=500GB
#SBATCH -c 8
#SBATCH --time=02-00:00:00
#Define sdout path
#SBATCH --output=/work/projects/age_sorter/scripts/GSE254569/output_GSE254569_2Seur.txt
#Define sderr path
#SBATCH --error=/work/projects/age_sorter/scripts/GSE254569/error_GSE254569_2Seur.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

Rscript h5ad_to_seur.R

conda deactivate