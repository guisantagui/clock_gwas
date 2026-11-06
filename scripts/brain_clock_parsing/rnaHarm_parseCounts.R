################################################################################
# Brain clock: Merge RNAseq harmonization counts files                         #
################################################################################

if (!require(plotUtils)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}

library(plotUtils)

# Directory stuff
################################################################################
setwd("../..")

outDir <- "./results/parsed/"

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

# Read RNAseqHarmonization and LBP files
################################################################################

RNAseqHarmDir <- "./data/expression/RNAseqHarm/"

rnaHarm_intMetDatFile <- "./results/metadata_parsed/RNAseq_Harmonization_ind_all.csv"

indMetDat <- read.csv(rnaHarm_intMetDatFile, row.names = 1)

ROSMAP_metDat <- "./data/metadata/ROSMAP/ROSMAP_clinical.csv"

ROSMAP_metDat <- read.csv(ROSMAP_metDat, row.names = 1)

Mayo <- read_table_fast(paste0(RNAseqHarmDir,
                               "Mayo_gene_all_counts_matrix_clean.txt"))

MSBB <- read_table_fast(paste0(RNAseqHarmDir,
                               "MSBB_gene_all_counts_matrix_clean.txt"))

ROSMAP1 <- read_table_fast(paste0(RNAseqHarmDir,
                                  "ROSMAP_batch1_gene_all_counts_matrix_clean.txt"))

ROSMAP2 <- read_table_fast(paste0(RNAseqHarmDir,
                                  "ROSMAP_batch2_gene_all_counts_matrix_clean.txt"))

ROSMAP3 <- read_table_fast(paste0(RNAseqHarmDir,
                                  "ROSMAP_batch3_gene_all_counts_matrix_clean.txt"))

ROSMAP4 <- read_table_fast(paste0(RNAseqHarmDir,
                                  "ROSMAP_batch4_gene_all_counts_matrix_clean.txt"))

# Inspect the datasets and compare to metadata
################################################################################

rosmapColnames <- c(colnames(ROSMAP1),
                    colnames(ROSMAP2),
                    colnames(ROSMAP3),
                    colnames(ROSMAP4))

all(make.names(indMetDat$specimenID[indMetDat$substudy == "ROSMAP"]) %in% rosmapColnames)
rosmapColnames[!rosmapColnames %in% make.names(indMetDat$specimenID)]

# There are some samples in the counts data that are not in the integrated
# metadata. Check individual metadata files to see if by any chance 
# These samples are the 90+



# All of the ROSMAP samples in integrated metadata are in ROSMAP counts files.

all(make.names(indMetDat$specimenID[indMetDat$substudy == "Mayo"]) %in% colnames(Mayo))

notInCounDF_mayo <- indMetDat$specimenID[indMetDat$substudy == "Mayo"][!make.names(indMetDat$specimenID[indMetDat$substudy == "Mayo"]) %in% colnames(Mayo)]
indMetDat[indMetDat$specimenID %in% notInCounDF_mayo, ]

# There are 21 samples from Mayo dataset that are not in the counts file. None of
# them are controls, so we migth just exclude them

all(make.names(indMetDat$specimenID[indMetDat$substudy == "MSBB"]) %in% colnames(MSBB))

notInCounDF_MSBB <- indMetDat$specimenID[indMetDat$substudy == "MSBB"][!make.names(indMetDat$specimenID[indMetDat$substudy == "MSBB"]) %in% colnames(MSBB)]
indMetDat[indMetDat$specimenID %in% notInCounDF_MSBB, ]

# There are 2 samples from MSBB dataset that are not in the counts file. None of
# them are controls, so we migth just exclude them

# The features in the 4 rosmap datasets are the same
all((ROSMAP1$feature == ROSMAP2$feature) == (ROSMAP3$feature == ROSMAP4$feature))

# Parse the datasets to have them in the same format
################################################################################

# Remove first four rows of the datasets, as they contain info about the
# mapping.

ROSMAP1 <- ROSMAP1[5:nrow(ROSMAP1), ]
ROSMAP2 <- ROSMAP2[5:nrow(ROSMAP2), ]
ROSMAP3 <- ROSMAP3[5:nrow(ROSMAP3), ]
ROSMAP4 <- ROSMAP4[5:nrow(ROSMAP4), ]

rownames(ROSMAP1) <- ROSMAP1$feature
rownames(ROSMAP2) <- ROSMAP2$feature
rownames(ROSMAP3) <- ROSMAP3$feature
rownames(ROSMAP4) <- ROSMAP4$feature

ROSMAP1 <- ROSMAP1[, colnames(ROSMAP1) != "feature"]
ROSMAP2 <- ROSMAP2[, colnames(ROSMAP2) != "feature"]
ROSMAP3 <- ROSMAP3[, colnames(ROSMAP3) != "feature"]
ROSMAP4 <- ROSMAP4[, colnames(ROSMAP4) != "feature"]

Mayo <- Mayo[5:nrow(Mayo), ]
rownames(Mayo) <- Mayo$feature
Mayo <- Mayo[, colnames(Mayo) != "feature"]

MSBB <- MSBB[5:nrow(MSBB), ]
rownames(MSBB) <- MSBB$feature
MSBB <- MSBB[, colnames(MSBB) != "feature"]

# Convert to TPMs
################################################################################

# Convert counts to TPMs
#all(rownames(Mayo) == rownames(MSBB))
#all(rownames(Mayo) == rownames(ROSMAP1))
#all(rownames(Mayo) == rownames(ROSMAP2))
#all(rownames(Mayo) == rownames(ROSMAP3))
#all(rownames(Mayo) == rownames(ROSMAP4))

# As the rownames (the features) of all the datasets is the same we only 
# need one geneInfo dataframe.

# Obtain gene info dataframe with bioMart
#human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# This is an alternative server for moments where the main might not work
#human <- useEnsembl(biomart = "genes", dataset="hsapiens_gene_ensembl", mirror = "asia")

#geneInfo <- getENSEMBL_ID_info(rownames(ROSMAP1), mart = human)


#geneInfo[geneInfo$ensembl_gene_id_version_new == geneInfo$ensembl_gene_id_version_new[duplicated(geneInfo$ensembl_gene_id_version_new)][3], ]
# Duplicated genes happen because there are copies in Y chromosome.

# Obtain TPMs of the datasets
#ROSMAP1_TPMs <- counts2TPMs(ROSMAP1, geneInfo = geneInfo)
#ROSMAP2_TPMs <- counts2TPMs(ROSMAP2, geneInfo = geneInfo)
#ROSMAP3_TPMs <- counts2TPMs(ROSMAP3, geneInfo = geneInfo)
#ROSMAP4_TPMs <- counts2TPMs(ROSMAP4, geneInfo = geneInfo)
#Mayo_TPMs <- counts2TPMs(Mayo, geneInfo = geneInfo)
#MSBB_TPMs <- counts2TPMs(MSBB, geneInfo = geneInfo)

# Add a column indicating ROSMAP batch to the integrated metadata.
################################################################################

indMetDat$ROSMAP_batch <- rep(NA, nrow(indMetDat))
#indMetDat$ROSMAP_batch[make.names(indMetDat$specimenID) %in% colnames(ROSMAP1_TPMs)] <- 1
#indMetDat$ROSMAP_batch[make.names(indMetDat$specimenID) %in% colnames(ROSMAP2_TPMs)] <- 2
#indMetDat$ROSMAP_batch[make.names(indMetDat$specimenID) %in% colnames(ROSMAP3_TPMs)] <- 3
#indMetDat$ROSMAP_batch[make.names(indMetDat$specimenID) %in% colnames(ROSMAP4_TPMs)] <- 4

indMetDat$ROSMAP_batch[make.names(indMetDat$specimenID) %in% colnames(ROSMAP1)] <- 1
indMetDat$ROSMAP_batch[make.names(indMetDat$specimenID) %in% colnames(ROSMAP2)] <- 2
indMetDat$ROSMAP_batch[make.names(indMetDat$specimenID) %in% colnames(ROSMAP3)] <- 3
indMetDat$ROSMAP_batch[make.names(indMetDat$specimenID) %in% colnames(ROSMAP4)] <- 4

intMetDat_wROSMAP_batchFile <- gsub("all", "all_ROSMAPBtch", rnaHarm_intMetDatFile)

indMetDat$mmse30_first_ad_dx <- rep(NA, nrow(indMetDat))
indMetDat$mmse30_lv <- rep(NA, nrow(indMetDat))

indMetDat$mmse30_first_ad_dx[indMetDat$substudy == "ROSMAP"] <- ROSMAP_metDat$cts_mmse30_first_ad_dx[match(indMetDat$individualID[indMetDat$substudy == "ROSMAP"],
                                                                                                           ROSMAP_metDat$individualID)]

indMetDat$mmse30_lv[indMetDat$substudy == "ROSMAP"] <- ROSMAP_metDat$cts_mmse30_lv[match(indMetDat$individualID[indMetDat$substudy == "ROSMAP"],
                                                                                         ROSMAP_metDat$individualID)]

write.csv(indMetDat, file = intMetDat_wROSMAP_batchFile)

# Do the transpose of the TPMs datasets and merge them.
#ROSMAP1_TPMs <- t(ROSMAP1_TPMs)
#ROSMAP2_TPMs <- t(ROSMAP2_TPMs)
#ROSMAP3_TPMs <- t(ROSMAP3_TPMs)
#ROSMAP4_TPMs <- t(ROSMAP4_TPMs)
#Mayo_TPMs <- t(Mayo_TPMs)
#MSBB_TPMs <- t(MSBB_TPMs)

# Merge the datasets in a single DF
#RNAseqHarm_allTPMs <- rbind(ROSMAP1_TPMs,
#                            ROSMAP2_TPMs,
#                            ROSMAP3_TPMs,
#                            ROSMAP4_TPMs,
#                            Mayo_TPMs,
#                            MSBB_TPMs)

# Merge counts datasets
################################################################################
RNAseqHarm_allCounts <- cbind(ROSMAP1,
                              ROSMAP2,
                              ROSMAP3,
                              ROSMAP4,
                              Mayo,
                              MSBB)

# Keep the samples that we already selected in the integrated RNAseq
# Harmonization metadata and write the result in a single file.

#RNAseqHarm_allTPMs <- RNAseqHarm_allTPMs[rownames(RNAseqHarm_allTPMs) %in% make.names(indMetDat$specimenID), ]

RNAseqHarm_allCounts <- RNAseqHarm_allCounts[, colnames(RNAseqHarm_allCounts) %in% make.names(indMetDat$specimenID)]

#write.csv(RNAseqHarm_allTPMs, paste0(outDir, "RNAseqHarm_allTPMs.csv"))
write.csv(RNAseqHarm_allCounts, paste0(outDir, "RNAseqHarm_allCounts.csv"))

plot(density(log(RNAseqHarm_allTPMs)))

#LBP_TPMs <- counts2TPMs(LBP,
#                        geneInfo = geneInfo_LBP,
#                        IDsIn = "ensembl_gene_id_version_orig")

#sum(gsub("\\..*", "", rownames(LBP)) %in% gsub("\\..*", "", colnames(RNAseqHarm_allCounts)))/nrow(LBP_TPMs)
# 99.9% of the genes in LBP_TPMs dataset are in the RNAseqHarm one.
#sum(gsub("\\..*", "", colnames(RNAseqHarm_allTPMs)) %in% gsub("\\..*", "", rownames(LBP_TPMs)))/ncol(RNAseqHarm_allTPMs)
# 97% of the genes in RNAseqHarm dataset are in the LBP_TPMs one.

# The overlap of RNAseqHarm and LBP is quite significant in terms of genome
# covered, so we can merge them.

# Go through LBP metadata to give it the shape we used for generating the
# integrated RNAseq Harmonization metadata

#LBP_parsedMetDat <- read.csv("./results/metadata_parsed/LBP_ind_all.csv")

#LBP_metDatInd <- read.csv("./data/metadata/LBP/LBP_individual_metadata.csv")
#LBP_metDatAss <- read.csv("./data/metadata/LBP/LBP_assay_RNAseq_metadata.csv")
#LBP_metDatBiosp <- read.csv("./data/metadata/LBP/LBP_biospecimen_metadata.csv")

#LBP_metDatInt <- data.frame(specimenID = LBP_metDatBiosp$specimenID[make.names(LBP_metDatBiosp$specimenID) %in% colnames(LBP_TPMs)])
#LBP_metDatInt$platform <- LBP_metDatAss$platform[match(LBP_metDatInt$specimenID,
#                                                       LBP_metDatAss$specimenID)]

#LBP_metDatInt$RIN <- LBP_metDatAss$RIN[match(LBP_metDatInt$specimenID,
#                                             LBP_metDatAss$specimenID)]

#LBP_metDatInt$libraryPrep <- LBP_metDatAss$libraryPrep[match(LBP_metDatInt$specimenID,
#                                                             LBP_metDatAss$specimenID)]

#LBP_metDatInt$libraryPreparationMethod <- LBP_metDatAss$libraryPreparationMethod[match(LBP_metDatInt$specimenID,
#                                                                                       LBP_metDatAss$specimenID)]

#LBP_metDatInt$runType <- LBP_metDatAss$runType[match(LBP_metDatInt$specimenID,
#                                                     LBP_metDatAss$specimenID)]

#LBP_metDatInt$readLength <- LBP_metDatAss$readLength[match(LBP_metDatInt$specimenID,
#                                                           LBP_metDatAss$specimenID)]

#LBP_metDatInt$individualID <- LBP_metDatBiosp$individualID[make.names(LBP_metDatBiosp$specimenID) %in% colnames(LBP_TPMs)]
#LBP_metDatInt$organ <- LBP_metDatBiosp$organ[make.names(LBP_metDatBiosp$specimenID) %in% colnames(LBP_TPMs)]
#LBP_metDatInt$tissue <- LBP_metDatBiosp$tissue[make.names(LBP_metDatBiosp$specimenID) %in% colnames(LBP_TPMs)]
#LBP_metDatInt$assay <- LBP_metDatBiosp$assay[make.names(LBP_metDatBiosp$specimenID) %in% colnames(LBP_TPMs)]
#LBP_metDatInt$exclude <- rep(NA, nrow(LBP_metDatInt))
#LBP_metDatInt$excludeReason <- rep("", nrow(LBP_metDatInt))
#LBP_metDatInt$sex <- LBP_metDatInd$sex[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
#LBP_metDatInt$race <- LBP_metDatInd$race[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
#LBP_metDatInt$ageDeath <- LBP_metDatInd$ageDeath[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
#LBP_metDatInt$apoeGenotype <- rep(NA, nrow(LBP_metDatInt))
#LBP_metDatInt$pmi <- LBP_metDatInd$pmi[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
#LBP_metDatInt$Braak <- LBP_metDatInd$Braak[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
#LBP_metDatInt$diagn_4BrainClck <- rep(NA, nrow(LBP_metDatInt))
#LBP_metDatInt$diagn_4BrainClck[match(LBP_parsedMetDat$individualID, LBP_metDatInt$individualID)]

#LBP_parsedMetDat$diagn_4BrainClck[match(LBP_metDatInt$individualID, LBP_parsedMetDat$individualID)]
#LBP_parsedMetDat$individualID[match(LBP_metDatInt$individualID, LBP_parsedMetDat$individualID)]

# Add the categories we built when parsing

#LBP_IDs_inParsedMetDat <- LBP_metDatInt$individualID[LBP_metDatInt$individualID %in% LBP_parsedMetDat$individualID]

#LBP_metDatInt$diagn_4BrainClck[LBP_metDatInt$individualID %in% LBP_parsedMetDat$individualID] <- LBP_parsedMetDat$diagn_4BrainClck[match(LBP_IDs_inParsedMetDat, LBP_parsedMetDat$individualID)]

#LBP_metDatInt$substudy <- rep("LBP", nrow(LBP_metDatInt))

#LBP_metDatInt$batch_rna <- LBP_metDatAss$rnaBatch[match(LBP_metDatInt$specimenID,
#                                                        LBP_metDatAss$specimenID)]

#LBP_metDatInt$batch_lib <- LBP_metDatAss$libraryBatch[match(LBP_metDatInt$specimenID,
#                                                            LBP_metDatAss$specimenID)]

#LBP_metDatInt$batch_seq <- LBP_metDatAss$sequencingBatch[match(LBP_metDatInt$specimenID,
#                                                               LBP_metDatAss$specimenID)]

#LBP_metDatInt <- LBP_metDatInt[!is.na(LBP_metDatInt$diagn_4BrainClck), ]
#table(LBP_metDatInt$ageDeath)

# Remove samples with age == ""
#LBP_metDatInt <- LBP_metDatInt[LBP_metDatInt$ageDeath != "", ]
#LBP_metDatInt$ageDeath <- as.numeric(LBP_metDatInt$ageDeath)

#write.csv(LBP_metDatInt, file = paste0(outDir, "LBP_metadata_unified.csv"))

#LBP_TPMs <- LBP_TPMs[, make.names(colnames(LBP_TPMs)) %in% LBP_metDatInt$specimenID]

write.csv(LBP_TPMs, file = paste0(outDir, "LBP_TPMs.csv"))

# Save geneInfo DFs
write.csv(geneInfo, paste0(outDir, "geneInfo_rnaSeqHarm.csv"))
write.csv(geneInfo_LBP, paste0(outDir, "geneInfo_LBP.csv"))

# Load GTEx data.

GTEx_TPMs <- read.csv(paste0(outDir, "GTEx_allTPMs.csv"), row.names = 1)

geneInfo_GTEx <- getENSEMBL_ID_info(GTEx_TPMs$Name,
                                    mart = human,
                                    filt = "ensembl_gene_id_version")

GTEx_TPMs[!GTEx_TPMs$Name %in% geneInfo_GTEx$ensembl_gene_id_version_orig, 1:3]

write.csv(geneInfo_GTEx, paste0(outDir, "geneInfo_GTEx.csv"))