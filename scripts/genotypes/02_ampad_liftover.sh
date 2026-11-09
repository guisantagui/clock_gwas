#!/bin/bash -l
################################################################################
# Genotypes, step 2: AMP-AD WGS (ROSMAP/Mayo/MSBB, b37) -> GRCh38.             #
#                                                                              #
# Per chromosome:                                                              #
#   1. keep PASS (VQSR) biallelic SNVs, set GQ<20 | DP<10 genotypes missing,  #
#      keep GT only (much smaller files for Picard)                            #
#   2. rename contigs 1..22,X -> chr1..chr22,chrX                              #
#   3. Picard LiftoverVcf with the UCSC hg19->hg38 chain onto the Broad        #
#      GRCh38 FASTA (REF/ALT swaps recovered; failures written to reject VCF) #
#   4. drop variants that moved to another chromosome, sort, then the common  #
#      standardisation (normalise vs FASTA, chr:pos:ref:alt IDs)              #
#                                                                              #
# Usage: sbatch --array=1-23 scripts/genotypes/02_ampad_liftover.sh            #
################################################################################
#SBATCH --job-name=ampad_liftover
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=32GB
#SBATCH --time=0-12:00:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH --output=logs/ampad_liftover_%a.out

set -euo pipefail
source "$(dirname "$0")/config.sh" 2>/dev/null || source scripts/genotypes/config.sh
activate_env

i=${SLURM_ARRAY_TASK_ID:?run as a SLURM array 1-23}
c=${CHROMS[$((i - 1))]}
in=$AMPAD_DIR/${AMPAD_PREFIX}_${c}.recalibrated_variants.vcf.gz
od=$OUT/ampad; tmp=$od/tmp_chr$c
mkdir -p "$od" "$tmp" logs

echo "[$(date)] chr$c: filter + rename"
for k in $(seq 1 22) X Y MT; do echo -e "$k\tchr$k"; done > "$tmp/rename.txt"
bcftools view -h "$in" > "$tmp/header.txt"
setgt="cat"
grep -q "ID=GQ," "$tmp/header.txt" && setgt="bcftools +setGT -Ou - -- -t q -n . -i FMT/GQ<${MIN_GQ}|FMT/DP<${MIN_DP}"
bcftools view -Ou -f PASS -m2 -M2 -v snps "$in" |
        $setgt |
        bcftools annotate -Ou -x INFO,^FORMAT/GT --rename-chrs "$tmp/rename.txt" - |
        bcftools view -Oz -o "$tmp/b37.vcf.gz" --threads "$THREADS" -
bcftools index -t "$tmp/b37.vcf.gz"

echo "[$(date)] chr$c: liftover"
picard -Xmx28g LiftoverVcf \
        I="$tmp/b37.vcf.gz" O="$tmp/hg38.vcf.gz" REJECT="$tmp/rejected.vcf.gz" \
        CHAIN="$CHAIN" R="$FASTA" \
        RECOVER_SWAPPED_REF_ALT=true WARN_ON_MISSING_CONTIG=true \
        MAX_RECORDS_IN_RAM=200000 TMP_DIR="$tmp"

echo "[$(date)] chr$c: same-chromosome, sort, standardise"
bcftools view -h "$tmp/hg38.vcf.gz" > "$tmp/header38.txt"
bcftools view -Ou -t "chr$c" "$tmp/hg38.vcf.gz" |
        bcftools sort -Ou -T "$tmp" - |
        standardise "$od/ampad_chr$c.bcf" "$tmp/header38.txt"

n_in=$(bcftools index -n "$tmp/b37.vcf.gz" 2>/dev/null || bcftools view -H "$tmp/b37.vcf.gz" | wc -l)
n_rej=$(bcftools view -H "$tmp/rejected.vcf.gz" | wc -l)
n_out=$(bcftools index -n "$od/ampad_chr$c.bcf")
echo -e "chr$c\t$n_in\t$n_rej\t$n_out" | tee "$od/liftover_stats_chr$c.tsv"
rm -rf "$tmp"
echo "[$(date)] chr$c done: in $n_in, rejected $n_rej, kept $n_out"
