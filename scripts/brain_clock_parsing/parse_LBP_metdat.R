################################################################################
# Brain clock: Unify LBP harmonization metadata file with RNAseq Harmonization #
# format                                                                       #
################################################################################

# Directory stuff
################################################################################
setwd("../..")

outDir <- "./results/parsed/"

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

# Load data
################################################################################
LBP_parsedMetDat <- read.csv("./results/metadata_parsed/LBP_ind_all.csv")

LBP_metDatInd <- read.csv("./data/metadata/LBP/LBP_individual_metadata.csv")
LBP_metDatAss <- read.csv("./data/metadata/LBP/LBP_assay_RNAseq_metadata.csv")
LBP_metDatBiosp <- read.csv("./data/metadata/LBP/LBP_biospecimen_metadata.csv")

# Build integrated metadata dataframe
################################################################################
LBP_metDatInt <- data.frame(specimenID = LBP_metDatBiosp$specimenID)

LBP_metDatInt$platform <- LBP_metDatAss$platform[match(LBP_metDatInt$specimenID,
                                                       LBP_metDatAss$specimenID)]

LBP_metDatInt$RIN <- LBP_metDatAss$RIN[match(LBP_metDatInt$specimenID,
                                             LBP_metDatAss$specimenID)]

LBP_metDatInt$libraryPrep <- LBP_metDatAss$libraryPrep[match(LBP_metDatInt$specimenID,
                                                             LBP_metDatAss$specimenID)]

LBP_metDatInt$libraryPreparationMethod <- LBP_metDatAss$libraryPreparationMethod[match(LBP_metDatInt$specimenID,
                                                                                       LBP_metDatAss$specimenID)]

LBP_metDatInt$runType <- LBP_metDatAss$runType[match(LBP_metDatInt$specimenID,
                                                     LBP_metDatAss$specimenID)]

LBP_metDatInt$readLength <- LBP_metDatAss$readLength[match(LBP_metDatInt$specimenID,
                                                           LBP_metDatAss$specimenID)]

LBP_metDatInt$individualID <- LBP_metDatBiosp$individualID
LBP_metDatInt$organ <- LBP_metDatBiosp$organ
LBP_metDatInt$tissue <- LBP_metDatBiosp$tissue
LBP_metDatInt$assay <- LBP_metDatBiosp$assay
LBP_metDatInt$exclude <- rep(NA, nrow(LBP_metDatInt))
LBP_metDatInt$excludeReason <- rep("", nrow(LBP_metDatInt))
LBP_metDatInt$sex <- LBP_metDatInd$sex[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
LBP_metDatInt$race <- LBP_metDatInd$race[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
LBP_metDatInt$ageDeath <- LBP_metDatInd$ageDeath[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
LBP_metDatInt$apoeGenotype <- rep(NA, nrow(LBP_metDatInt))
LBP_metDatInt$pmi <- LBP_metDatInd$pmi[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
LBP_metDatInt$Braak <- LBP_metDatInd$Braak[match(LBP_metDatInt$individualID, LBP_metDatInd$individualID)]
LBP_metDatInt$diagn_4BrainClck <- rep(NA, nrow(LBP_metDatInt))
LBP_metDatInt$diagn_4BrainClck[match(LBP_parsedMetDat$individualID, LBP_metDatInt$individualID)]

LBP_parsedMetDat$diagn_4BrainClck[match(LBP_metDatInt$individualID, LBP_parsedMetDat$individualID)]
LBP_parsedMetDat$individualID[match(LBP_metDatInt$individualID, LBP_parsedMetDat$individualID)]

# Add the categories we built when parsing

LBP_metDatInt <- LBP_metDatInt[LBP_metDatInt$individualID %in% LBP_parsedMetDat$individualID, ]


LBP_metDatInt$diagn_4BrainClck <- LBP_parsedMetDat$diagn_4BrainClck[match(LBP_metDatInt$individualID, LBP_parsedMetDat$individualID)]

LBP_metDatInt$substudy <- rep("LBP", nrow(LBP_metDatInt))

LBP_metDatInt$batch_rna <- LBP_metDatAss$rnaBatch[match(LBP_metDatInt$specimenID,
                                                        LBP_metDatAss$specimenID)]

LBP_metDatInt$batch_lib <- LBP_metDatAss$libraryBatch[match(LBP_metDatInt$specimenID,
                                                            LBP_metDatAss$specimenID)]

LBP_metDatInt$batch_seq <- LBP_metDatAss$sequencingBatch[match(LBP_metDatInt$specimenID,
                                                               LBP_metDatAss$specimenID)]

LBP_metDatInt <- LBP_metDatInt[!is.na(LBP_metDatInt$diagn_4BrainClck), ]
table(LBP_metDatInt$ageDeath)

# Remove samples with age == ""
LBP_metDatInt <- LBP_metDatInt[LBP_metDatInt$ageDeath != "", ]
LBP_metDatInt$ageDeath <- as.numeric(LBP_metDatInt$ageDeath)

write.csv(LBP_metDatInt, file = paste0(outDir, "LBP_metadata_unified.csv"))