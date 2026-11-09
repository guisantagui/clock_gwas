#!/bin/bash -l
################################################################################
# Genotypes, step 0: inventory of the WGS callsets before merging.            #
# Prints build, sample counts/IDs, file layout and contents of each callset.  #
#                                                                              #
# Usage (on a compute node, gwas env active):                                  #
#   bash scripts/genotypes/00_inventory.sh > inventory.txt 2>&1                #
################################################################################

G=${G:-/scratch/users/gsantamaria/test_large_files/gwas_clock}
AMPAD=$G/AMP_AD/WGS_Harmonization/genomics
LBP=$G/AMP_AD/LBP/alignments
GTEX=$G/GTEx/genomics
NABEC=$G/NABEC/genomics

echo "### Tools"
bcftools --version | head -1; plink2 --version; gatk --version 2>/dev/null | head -1

echo; echo "### AMP-AD (ROSMAP/Mayo/MSBB) joint VCF"
f=$(ls $AMPAD/*_22.recalibrated_variants.vcf.gz)
bcftools view -h "$f" | grep -m3 "^##contig"
bcftools view -h "$f" | grep -m1 "^##reference"
echo "n samples: $(bcftools query -l "$f" | wc -l)"
echo "first sample IDs:"; bcftools query -l "$f" | head -5
echo "FILTER values (first 200k records chr22):"
bcftools query -f '%FILTER\n' "$f" | head -200000 | sort | uniq -c | sort -rn | head
echo "first records:"; bcftools view -H "$f" | head -2 | cut -c1-300

echo; echo "### LBP Sarek output"
echo "n sample dirs: $(ls $LBP/variant_calling/haplotypecaller | wc -l)"
ls -la $LBP/variant_calling/haplotypecaller/5721-EC-1/
echo "samples missing a gVCF:"
for d in $LBP/variant_calling/haplotypecaller/*/; do
        s=$(basename "$d"); ls "$d"/*.g.vcf.gz >/dev/null 2>&1 || echo "  $s"
done
gv=$(ls $LBP/variant_calling/haplotypecaller/5721-EC-1/*.g.vcf.gz 2>/dev/null | head -1)
[[ -n "$gv" ]] && { bcftools view -h "$gv" | grep -m2 "^##contig"; bcftools query -l "$gv"; }
echo "preprocessing:"; ls $LBP/preprocessing/ 2>/dev/null; ls $LBP/preprocessing/* 2>/dev/null | head -5
echo "sarek reference:"; ls $LBP/reference/ 2>/dev/null | head; ls $G/AMP_AD/LBP/reference/ 2>/dev/null | head

echo; echo "### GTEx"
ls -la $GTEX
tar -tvf $GTEX/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1.GRU.tar | head -20

echo; echo "### NABEC"
for t in $NABEC/phg00*.tar.gz; do echo "== $t"; tar -tzvf "$t" | head -20; done
ls $NABEC/phs000249 $NABEC/phs001462 2>/dev/null | head -20

echo; echo "### Disk"
df -h $G | tail -1
