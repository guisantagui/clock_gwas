################################################################################
# Brain clock: combine the integrated dataset made up with RNAseqHarm, LBP,    #
# GTEx and TBI with brainSeq phase II, AgeAnno pseudo-bulk data and our        #
# perturbation compilation.                                                    #
################################################################################

if (!require("devtools",quietly = T)){
    install.packages("devtools",
                     repos = 'http://cran.us.r-project.org')
}

if(!require("plotUtils", quietly = T)){
        devtools::install_github('guisantagui/plotUtils', upgrade = "never")
}
library(plotUtils)
library(dplyr)
# Directory stuff
################################################################################

comb_counts_f <- "../../results/parsed/combined_counts_wTBI.csv"
comb_metdat_f <- "../../results/parsed/combined_metDat_wTBI.csv"

bs_pII_counts_f <- "../../results/parsed/brainseq_pII_counts.csv"
bs_pII_metdat_f <- "../../results/parsed/brainseq_pII_metadata.csv"

bs_pI_counts_f <- "../../results/parsed/brainseq_pI_counts.csv"
bs_pI_metdat_f <- "../../results/parsed/brainseq_pI_metadata.csv"

aa_counts_f <- "../../results/parsed/ageAnno/ageanno_brain_pb_counts.csv"
aa_metdat_f <- "../../results/parsed/ageAnno/ageanno_brain_pb_metdat.csv"

#pert_counts_f <- "./results/parsed/counts_111.csv"
#pert_metdat_f <- "./results/parsed/metDat_111.csv"

out_dir <- "../../results/parsed/merged/"
create_dir_if_not(out_dir)

# Load data
################################################################################
comb_counts <- read_table_fast(comb_counts_f, row.names = 1)
comb_metdat <- read_table_fast(comb_metdat_f, row.names = 1)

bs_pII_counts <- read_table_fast(bs_pII_counts_f, row.names = 1)
bs_pII_metdat <- read_table_fast(bs_pII_metdat_f, row.names = 1)

bs_pI_counts <- read_table_fast(bs_pI_counts_f, row.names = 1)
bs_pI_metdat <- read_table_fast(bs_pI_metdat_f, row.names = 1)

aa_counts <- read_table_fast(aa_counts_f, row.names = 1)
aa_metdat <- read_table_fast(aa_metdat_f, row.names = 1)

#pert_counts <- read_table_fast(pert_counts_f, row.names = 1)
#pert_metdat <- read_table_fast(pert_metdat_f, row.names = 1)

# Merge counts
################################################################################

common_genes <- colnames(comb_counts) %>%
        intersect(rownames(bs_pII_counts)) %>%
        intersect(rownames(bs_pI_counts)) %>%
        intersect(rownames(aa_counts))# %>%
        #intersect(rownames(pert_counts))

bs_pII_counts <- t(bs_pII_counts)
bs_pI_counts <- t(bs_pI_counts)
aa_counts <- t(aa_counts)
#pert_counts <- t(pert_counts)

merged_counts <- rbind.data.frame(
        comb_counts[, common_genes],
        bs_pII_counts[, common_genes],
        bs_pI_counts[, common_genes],
        aa_counts[, common_genes]
)

# Merge metadata
################################################################################

# Add perturbation columns
comb_metdat$perturbation <- "none"
comb_metdat$exper_group <- ""

cols2Add_bs_pII <- colnames(comb_metdat)[!colnames(comb_metdat) %in% colnames(bs_pII_metdat)]
cols2Add_bs_pI <- colnames(comb_metdat)[!colnames(comb_metdat) %in% colnames(bs_pI_metdat)]
cols2Add_aa <- colnames(comb_metdat)[!colnames(comb_metdat) %in% colnames(aa_metdat)]

bs_pII_metdat <- cbind.data.frame(bs_pII_metdat,
                                  data.frame(matrix(nrow = nrow(bs_pII_metdat),
                                                    ncol = length(cols2Add_bs_pII),
                                                    dimnames = list(rownames(bs_pII_metdat),
                                                                             cols2Add_bs_pII))))

bs_pI_metdat <- cbind.data.frame(bs_pI_metdat,
                                 data.frame(matrix(nrow = nrow(bs_pI_metdat),
                                                   ncol = length(cols2Add_bs_pI),
                                                   dimnames = list(rownames(bs_pI_metdat),
                                                                            cols2Add_bs_pI))))

aa_metdat <- cbind.data.frame(aa_metdat,
                              data.frame(matrix(nrow = nrow(aa_metdat),
                                                ncol = length(cols2Add_aa),
                                                dimnames = list(rownames(aa_metdat),
                                                                         cols2Add_aa))))

#pert_metdat <- cbind.data.frame(pert_metdat,
#                                data.frame(matrix(nrow = nrow(pert_metdat),
#                                                  ncol = length(cols2Add_pert),
#                                                  dimnames = list(rownames(pert_metdat),
#                                                                           cols2Add_pert))))

#pert_metdat$sex <- "undefined"

merged_metdat <- rbind.data.frame(
        comb_metdat,
        bs_pII_metdat[, colnames(comb_metdat)],
        bs_pI_metdat[, colnames(comb_metdat)],
        aa_metdat[, colnames(comb_metdat)]#,
        #pert_metdat[, colnames(comb_metdat)]
)

# Round ageDeath, as brainSeq has decimals
merged_metdat$ageDeath <- round(merged_metdat$ageDeath)

# Ensure that perturbation is set to "none" and exper group to ""
# in all the samples in the merged dataset (all are clinical samples)
merged_metdat$perturbation <- "none"
merged_metdat$exper_group <- ""

# Filter counts to keep only what is available in the metadata
merged_counts <- merged_counts[match(make.names(merged_metdat$specimenID),
                                     make.names(rownames(merged_counts))), ]

# Write the files
################################################################################

write_table_fast(merged_counts, sprintf("%smerged_counts.csv", out_dir))
write_table_fast(merged_metdat, sprintf("%smerged_metdat.csv", out_dir))