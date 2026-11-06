#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=dl_GSE254569_fastqs
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=02-00:00:00
#Define sdout path
#SBATCH --output=/work/projects/age_sorter/scripts/GSE254569/output_PA_rnaseq.txt
#Define sderr path
#SBATCH --error=/work/projects/age_sorter/scripts/GSE254569/error_PA_rnaseq.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate bioinfo
# Variables
####################################################################################################

accessions_file="/work/projects/age_sorter/data/GSE254569_sra_accessions/accessions.txt"
split_fastq=1
fastq_dir="/scratch/users/gsantamaria/test_larger_files/aging_trajs_data/fastq/GSE254569/"

mkdir -p "$fastq_dir"

cd $fastq_dir

download_accession() {
    accession=$1
    fastq_1="${accession}_1.fastq"
    fastq_2="${accession}_2.fastq"
    if [[ ! -f "$fastq_1" || ! -f "$fastq_2" ]]; then
        echo "Downloading accession: $accession"
        prefetch $accession
        fasterq-dump --threads 2 --split-files --include-technical $accession
    else
        echo "$accession has already been dumped. Skipping."
    fi
}

export -f download_accession
cat "$accessions_file" | parallel -j 4 download_accession {}

#while read -r accession; do
#  fastq_1="${accession}_1.fastq"
#  fastq_2="${accession}_2.fastq"
#  if [[ ! -f "$fastq_1" || ! -f "$fastq_2" ]]; then
#    echo "Downloading accession: $accession"
#    prefetch $accession
#    fasterq-dump --threads 2 --split-files --include-technical $accession
#  else
#    echo "$accession has already been dumped. Skipping."
#  fi
#done < $accessions_file

conda deactivate