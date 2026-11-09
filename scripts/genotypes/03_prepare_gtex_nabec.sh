#!/bin/bash -l
################################################################################
# Genotypes, step 3: GTEx v9 and NABEC WGS (both GRCh38) -> standardised BCF.  #
#                                                                              #
# GTEx: the 680 GB VCF is streamed straight out of the dbGaP tar (no need to  #
#       extract it), filtered to PASS biallelic SNVs with low-quality          #
#       genotypes set missing, GT only, and split per chromosome.             #
# NABEC: per-chromosome VCFs are extracted from the tar and standardised.     #
#                                                                              #
# Usage:                                                                       #
#   sbatch scripts/genotypes/03_prepare_gtex_nabec.sh gtex                     #
#   sbatch --array=1-23 scripts/genotypes/03_prepare_gtex_nabec.sh nabec       #
################################################################################
#SBATCH --job-name=prep_gtex_nabec
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=16GB
#SBATCH --time=2-00:00:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH --output=logs/prep_%x_%a.out

set -euo pipefail
source "$(dirname "$0")/config.sh" 2>/dev/null || source scripts/genotypes/config.sh
activate_env
mkdir -p logs
cohort=${1:?gtex or nabec}

if [[ $cohort == gtex ]]; then
        od=$OUT/gtex; mkdir -p "$od"
        echo "[$(date)] GTEx: header"
        tar -xOf "$GTEX_TAR" "$GTEX_VCF_IN_TAR" | bcftools view -h > "$od/header.txt" || true
        echo "[$(date)] GTEx: streaming filter (one pass over 680 GB)"
        tar -xOf "$GTEX_TAR" "$GTEX_VCF_IN_TAR" |
                bcftools view -Ou -f PASS,. -m2 -M2 -v snps - |
                standardise "$od/gtex_all.bcf" "$od/header.txt"
        echo "[$(date)] GTEx: split per chromosome"
        for c in "${CHROMS[@]}"; do
                bcftools view -Ob -o "$od/gtex_chr$c.bcf" -r "chr$c" "$od/gtex_all.bcf"
                bcftools index -f "$od/gtex_chr$c.bcf"
        done
        bcftools query -l "$od/gtex_all.bcf" > "$od/samples.txt"
        echo "[$(date)] GTEx done: $(wc -l < "$od/samples.txt") samples"

elif [[ $cohort == nabec ]]; then
        i=${SLURM_ARRAY_TASK_ID:?run nabec as a SLURM array 1-23}
        c=${CHROMS[$((i - 1))]}
        od=$OUT/nabec; tmp=$od/tmp_chr$c; mkdir -p "$od" "$tmp"
        f=$NABEC_DIR_IN_TAR/nabec.freeze9.chr$c.gtonly.vcf.gz
        tar -xzf "$NABEC_TAR" -C "$tmp" "$f"
        bcftools view -h "$tmp/$f" > "$tmp/header.txt"
        bcftools view -Ou -f PASS,. "$tmp/$f" | standardise "$od/nabec_chr$c.bcf" "$tmp/header.txt"
        [[ $c == 22 ]] && bcftools query -l "$od/nabec_chr$c.bcf" > "$od/samples.txt"
        rm -rf "$tmp"
        echo "[$(date)] NABEC chr$c done: $(bcftools index -n "$od/nabec_chr$c.bcf") variants"
else
        echo "unknown cohort $cohort"; exit 1
fi
