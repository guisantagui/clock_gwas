################################################################################
# Shared paths and settings for the genotype scripts. Sourced by every step.  #
################################################################################

G=/scratch/users/gsantamaria/test_large_files/gwas_clock
OUT=$G/genotypes                                   # all derived genotype files
REF=$G/reference
FASTA=$REF/Homo_sapiens_assembly38.fasta
CHAIN=$REF/hg19ToHg38.over.chain.gz

AMPAD_DIR=$G/AMP_AD/WGS_Harmonization/genomics
AMPAD_PREFIX=NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08
GTEX_TAR=$G/GTEx/genomics/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1.GRU.tar
GTEX_VCF_IN_TAR=./phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz
NABEC_TAR=$G/NABEC/genomics/phg001711.v1.NABEC_WGS.genotype-calls-vcf.WGS_markerset_grc38.c1.GRU.tar.gz
NABEC_DIR_IN_TAR=dbGaP_phs2636_NABEC_WGS_cohort_genotype_per_chromosome_VCF

# Chromosomes analysed (SLURM array index 1..23 -> chr1..chr22, chrX)
CHROMS=($(seq 1 22) X)

# Genotype-level filter: calls with low quality are set to missing.
MIN_GQ=20
MIN_DP=10

# Variant-level filters applied per cohort (before merging).
MAX_MISSING=0.05     # plink2 --geno
MIN_MAF=0.01         # plink2 --maf (GWAS of common variants)
HWE_P=1e-10          # plink2 --hwe (lenient: cohorts include cases)

THREADS=${SLURM_CPUS_PER_TASK:-4}

activate_env() {
        eval "$(micromamba shell hook --shell bash)"
        micromamba activate gwas
}

# Common tail of every per-cohort pipeline: keep biallelic SNVs, left-normalise
# against GRCh38 (records whose REF does not match the FASTA are dropped), set
# low-quality genotypes to missing (only if GQ/DP exist), keep GT only and give
# every variant the ID chr:pos:ref:alt. Reads VCF/BCF on stdin, writes BCF.
standardise() {
        local out=$1 hdr=$2
        local setgt="cat"
        if grep -q "ID=GQ," "$hdr" && grep -q "ID=DP," "$hdr"; then
                setgt="bcftools +setGT -Ou - -- -t q -n . -i FMT/GQ<${MIN_GQ}|FMT/DP<${MIN_DP}"
        fi
        bcftools view -Ou -m2 -M2 -v snps - |
                bcftools norm -Ou -f "$FASTA" -c x - |
                $setgt |
                bcftools annotate -Ou -x INFO,^FORMAT/GT --set-id '%CHROM:%POS:%REF:%ALT' - |
                bcftools view -Ob -o "$out" --threads "$THREADS" -
        bcftools index -f "$out"
}
