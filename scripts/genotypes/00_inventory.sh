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

command -v bcftools >/dev/null || { echo "bcftools not found: run 'micromamba activate gwas' first"; exit 1; }

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
HC=$LBP/variant_calling/haplotypecaller
echo "n sample dirs: $(ls $HC | wc -l)"
v=$HC/5721-EC-1/5721-EC-1.haplotypecaller.vcf.gz
echo "Is it a gVCF? (GVCFBlock header lines / <NON_REF> alleles):"
bcftools view -h "$v" | grep -c "^##GVCFBlock"
bcftools view -H "$v" | head -1000 | grep -c "NON_REF"
bcftools view -h "$v" | grep -m2 "^##contig"
bcftools view -h "$v" | grep -m1 "^##GATKCommandLine" | cut -c1-400
echo "n records: $(bcftools index -n $v)"; bcftools view -H "$v" | head -3 | cut -c1-250
echo "samples without a haplotypecaller VCF:"
for d in $HC/*/; do s=$(basename "$d"); [[ -s "$d/$s.haplotypecaller.vcf.gz" ]] || echo "  $s"; done
echo "CRAM/BAM after markduplicates:"
MD=$LBP/preprocessing/markduplicates
echo "n dirs: $(ls $MD | wc -l)"; ls -la $MD/5721-EC-1/ 2>/dev/null || ls -la $MD/$(ls $MD | head -1)/
echo "n CRAMs: $(ls $MD/*/*.cram 2>/dev/null | wc -l)"
echo "reference used by Sarek:"; ls -R $LBP/reference/ | head -20
find /scratch/users/gsantamaria/test_large_files/gwas_clock/AMP_AD/LBP -maxdepth 3 -name "*.fasta*" -o -maxdepth 3 -name "*.fa" | head

echo; echo "### GTEx"
ls -la $GTEX
echo "(single 680 GB vcf.gz inside the tar; needs extraction)"

echo; echo "### NABEC"
tar -xzOf $NABEC/phg001711.v1.NABEC_WGS.genotype-calls-vcf.WGS_markerset_grc38.c1.GRU.tar.gz dbGaP_phs2636_NABEC_WGS_cohort_genotype_per_chromosome_VCF/nabec.freeze9.chr22.gtonly.vcf.gz | bcftools view -h | grep -E -m4 "^##contig|^##reference|^#CHROM" | cut -c1-300
tar -tzf $NABEC/phg001711.v1.NABEC_WGS.genotype-calls-vcf.WGS_markerset_grc38.c1.GRU.tar.gz | grep -v "chr[0-9]*\.gtonly" 

echo; echo "### Disk"
df -h $G | tail -1
