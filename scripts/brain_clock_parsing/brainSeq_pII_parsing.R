################################################################################
# Brain clock: parse brainSeq Phase II data (from                              #
# doi: 10.1016/j.neuron.2019.05.013)                                           #
################################################################################
setwd("../..")

gene_file <- "./data/brainSeq_pII/rse_gene_unfiltered.Rdata"
out_dir <- "./results/parsed/"

# Load the data
################################################################################
load(gene_file)

# Parse metadata
################################################################################
met_dat <- rse_gene@colData

# Keep only controls (non-controls are schizophrenia patients)
met_dat <- met_dat[met_dat$Dx == "Control", ]
dim(met_dat)

# Keep samples with age higher than 20 yo
met_dat <- met_dat[met_dat$Age >= 20, ]

# Create unified metadata dataframe, with same column names as the rest
met_dat_unified <- met_dat[, c("RNum",
                               "RIN",
                               "BrNum",
                               "Region",
                               "Sex",
                               "Race",
                               "Age",
                               "Dx")]
met_dat_unified$RIN <- unlist(met_dat_unified$RIN)
met_dat_unified <- as.data.frame(met_dat_unified)
met_dat_unified$platform <- "TruSeq.v1"
met_dat_unified$libraryPrep <- ""
met_dat_unified$libraryPreparationMethod <- ""
met_dat_unified$runType <- ""
met_dat_unified$readLength <- NA
met_dat_unified$organ <- "brain"
met_dat_unified$assay <- "rnaSeq"
met_dat_unified$exclude <- NA
met_dat_unified$excludeReason <- ""
met_dat_unified$Sex <- gsub("F", "female", met_dat_unified$Sex)
met_dat_unified$Sex <- gsub("M", "male", met_dat_unified$Sex)
met_dat_unified$Race <- gsub("AA", "Black", met_dat_unified$Race)
met_dat_unified$Race <- gsub("AS", "Asian", met_dat_unified$Race)
met_dat_unified$Race <- gsub("CAUC", "White", met_dat_unified$Race)
met_dat_unified$Race <- gsub("HISP", "Hispanic", met_dat_unified$Race)
met_dat_unified$apoeGenotype <- ""
met_dat_unified$pmi <- NA
met_dat_unified$Braak <- NA
met_dat_unified$substudy <- "brainSeq_pII"
met_dat_unified$batch_rna <- NA
met_dat_unified$batch_lib <- ""
met_dat_unified$batch_seq <- NA
met_dat_unified$mmse30_first_ad_dx <- NA
met_dat_unified$mmse30_lv <- NA
met_dat_unified$perturbation <- "none"
met_dat_unified$exper_group <- ""

colnames(met_dat_unified) <- gsub("RNum", "specimenID", colnames(met_dat_unified))
colnames(met_dat_unified) <- gsub("BrNum", "individualID", colnames(met_dat_unified))
colnames(met_dat_unified) <- gsub("Region", "tissue", colnames(met_dat_unified))
colnames(met_dat_unified) <- gsub("Sex", "sex", colnames(met_dat_unified))
colnames(met_dat_unified) <- gsub("Race", "race", colnames(met_dat_unified))
colnames(met_dat_unified) <- gsub("Age", "ageDeath", colnames(met_dat_unified))
colnames(met_dat_unified) <- gsub("Dx", "diagn_4BrainClck", colnames(met_dat_unified))
met_dat_unified <- met_dat_unified[, c("specimenID",
                                       "platform",
                                       "RIN",
                                       "libraryPrep",
                                       "libraryPreparationMethod",
                                       "runType",
                                       "readLength",
                                       "individualID",
                                       "organ",
                                       "tissue",
                                       "assay",
                                       "exclude",
                                       "excludeReason",
                                       "sex",
                                       "race",
                                       "ageDeath",
                                       "apoeGenotype",
                                       "pmi",
                                       "Braak",
                                       "diagn_4BrainClck",
                                       "substudy",
                                       "batch_rna",
                                       "batch_lib",
                                       "batch_seq",
                                       "mmse30_first_ad_dx",
                                       "mmse30_lv",
                                       "perturbation",
                                       "exper_group")]
rownames(met_dat_unified) <- 1:nrow(met_dat_unified)

# Parse counts
################################################################################
counts <- rse_gene@assays$data$counts

# Change row names to ENSEMBL IDs using the rowRanges information
gene_info <- as.data.frame(rse_gene@rowRanges)
rownames(counts) <- gene_info$ensemblID[match(rownames(counts),
                                              gene_info$gencodeID)]
counts <- counts[, colnames(counts) %in% met_dat_unified$specimenID]
dup_genes <- rownames(counts)[duplicated(rownames(counts))]


# After changing to ENSEMBL IDs there are a few duplicated genes, so sum their
# counts
sum_dups <- function(m){
        dup_genes <- rownames(m)[duplicated(rownames(m))]
        for(i in seq_along(dup_genes)){
                dg <- dup_genes[i]
                nu_vec <- colSums(m[rownames(m) == dg, ])
                matches <- which(rownames(m) == dg)
                m[matches[1], ] <- nu_vec
                m <- m[-matches[-1], ]
        }
        return(m)
}

counts <- sum_dups(counts)

write.csv(counts, file = sprintf("%sbrainseq_pII_counts.csv",
                                 out_dir))
write.csv(met_dat_unified, file = sprintf("%sbrainseq_pII_metadata.csv",
                                          out_dir))