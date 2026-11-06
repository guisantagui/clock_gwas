

conda activate r-4.3.1

# 1. Parsing of AMP-AD datasets
# 1.1. Parsing of AMP-AD individual metadata
Rscript amp_ad_metadata_parsing.R

# 1.2. Homogenizing RNAseq Hamronization metadata and counts matrices
Rscript rnaHarm_parseCounts.R

# 1.3. Homogenizing LBP metadata
Rscript parse_LBP_metdat.R

# 2. Parsing of GTEx datasets
Rscript GTEx_parsing.R

# 3. Merging GTEx, LBP and RNAseq Harmonization datasets
Rscript merge_GTEx_LBP_rnaSeqHarm.R

# 4. Obtention of TBI data and merging with the combined dataset
Rscript tbi_get.R

# 5. Obtention of BrainSeq Phase 2 data
Rscript brainSeq_pII_parsing.R

# 6. Obtention of BrainSeq Phase 1 data
Rscript brainSeq_pI_parsing.R

# 7. Pre-processing, integration and annotation of AgeAnno brain scRNA-seq data
Rscript ageAnno_preproc.R

# 8. Integrate all clinical samples into a single file and pre-process
# 8.1 Add brainSeq Phase II, AgeAnno and perturbation compilation to merged dataset
Rscript merge_clinical.R

# 8.2 Filter clinical dataset
Rscript filt_dataset.R