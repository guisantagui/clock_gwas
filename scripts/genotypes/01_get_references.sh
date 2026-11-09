#!/bin/bash -l
################################################################################
# Genotypes, step 1: reference files.                                          #
#                                                                              #
# * GRCh38 FASTA used by Sarek (--genome GATK.GRCh38 = Broad                   #
#   Homo_sapiens_assembly38.fasta). Needed to decode the LBP CRAMs, to         #
#   left-normalise every callset against the same reference, and as the       #
#   liftover target.                                                           #
# * UCSC hg19 -> hg38 chain, for lifting the AMP-AD callset (b37, contigs     #
#   named 1..22,X) to GRCh38 after renaming its contigs to chr1..chr22,chrX.  #
#                                                                              #
# Usage: bash scripts/genotypes/01_get_references.sh [refdir]                  #
################################################################################
set -euo pipefail

REF=${1:-/scratch/users/gsantamaria/test_large_files/gwas_clock/reference}
mkdir -p "$REF"; cd "$REF"

BROAD=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
for f in Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict; do
        [[ -s $f ]] || wget -q --show-progress "$BROAD/$f"
done

[[ -s hg19ToHg38.over.chain.gz ]] || wget -q --show-progress \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# b37 contig names (1, 2, ..., X) -> UCSC style (chr1, ..., chrX)
for c in $(seq 1 22) X Y; do echo -e "$c\tchr$c"; done > b37_to_chr.txt

# Picard LiftoverVcf needs the target reference with a .dict next to it (done above).
md5sum Homo_sapiens_assembly38.fasta > Homo_sapiens_assembly38.fasta.md5
ls -la "$REF"
