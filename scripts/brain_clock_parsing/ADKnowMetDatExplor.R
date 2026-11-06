# Inspect metadata from AD knowledge portal studies

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/metadata_parsed/"

# AMP-AD_DiverseCohorts
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn51732482



if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}

AMP_AD_DiverseCohorts <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/AMP-AD_DiverseCohorts"

AMP_AD_DiverseCohorts_files <- list.files(AMP_AD_DiverseCohorts, full.names = T)

# Biospecimen
bspcmn <- read.csv(AMP_AD_DiverseCohorts_files[1])

# Assay
assay <- read.csv(AMP_AD_DiverseCohorts_files[2])

# Individual
indiv <- read.csv(AMP_AD_DiverseCohorts_files[3])

sum(indiv$ADoutcome == "Control")

grep("mayo", tolower(colnames(bspcmn)))
grep("mayo", tolower(colnames(assay)))
indiv[, grep("mayo", tolower(colnames(indiv)))]

# In indiv there is ADoutcome and mayoDx which assess if they are or not
# controls. Controls were assessed as not having AD based on braak <= Stage III
# and amyThal <= phase 2. ADOutcome variable is harmonized, so probably this is
# the good one. See dictionary for clarification. Most of the control samples
# according to this variable have ε3/ε3 genotype, with is average risk.
table(indiv[indiv$ADoutcome == "Control" & indiv$amyCerad == "None/No AD/C0", ]$apoeGenotype)

indiv[indiv$amyCerad == "None/No AD/C0", ]

indiv[indiv$ADoutcome == "Control" & indiv$amyCerad == "None/No AD/C0", ]

isCtrl <- indiv$ADoutcome == "Control" & indiv$amyCerad == "None/No AD/C0" & indiv$amyThal == "None"

indiv$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrl) + 1]



# Remove the Mayo, MSBB and ROSMAP cohorts, as we have them already covered
# in the RNAseq Harmonization study

table(indiv$cohort)

indiv <- indiv[!indiv$cohort %in% c("Mt Sinai Brain Bank", "ROS", "MAP", "Mayo Clinic"), ]

# Remove ambiguous age categories
table(indiv$ageDeath)
indiv <- indiv[indiv$ageDeath != "90+" & indiv$ageDeath != "missing or unknown", ]
indiv$ageDeath <- round(as.numeric(indiv$ageDeath))


indivCtrls <- indiv[indiv$diagn_4BrainClck == "Control", ]

dim(indivCtrls)

table(indiv[indiv$ADoutcome == "Control" & indiv$amyCerad == "None/No AD/C0" , ]$amyThal)

table(indivCtrls$ageDeath)

write.csv(indiv, paste0(outDir, "AMP_AD_ind_all.csv"))
write.csv(indivCtrls, paste0(outDir, "AMP_AD_ind_ctrls.csv"))

# rnaSeqReprocessing
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn9702085

rnaSeqReprocessing <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/rnaSeqReprocessing"
rnaSeqReprocessing_files <- list.files(rnaSeqReprocessing, full.names = T)

mayoIndiv <- read.csv(rnaSeqReprocessing_files[1])

mayoIndivCtrls <- mayoIndiv[mayoIndiv$diagnosis == "control" & mayoIndiv$Thal == 0, ]

mayoIndivCtrls <- mayoIndivCtrls[!is.na(mayoIndivCtrls$diagnosis), ]

mayoIndivCtrls <- mayoIndivCtrls[!is.na(mayoIndivCtrls$individualID), ]

# There is a category of controls and the category of Thal, which probably is
# amyloid in thalamus. Set to zero.

MSBBIndiv <- read.csv(rnaSeqReprocessing_files[2])
MSBBIndivCtrl <- MSBBIndiv[MSBBIndiv$CERAD == 4 & MSBBIndiv$Braak <= 3, ]
dim(MSBBIndivCtrl)

# No diagnosis values in MSBB. Use CERAD and Braak isntead

ROSMAPIndiv <- read.csv(rnaSeqReprocessing_files[4])


# dcfdx_lv is clinical diagnosis of cognitive status. 1 means NCI (no cognitive
# impairment). cogdx is the postmortem summary diagnosis. CERAD score is a
# semiquantitative metric measure of neuritic plates. Braak is something similar
ROSMAPIndivCtrls <- ROSMAPIndiv[ROSMAPIndiv$dcfdx_lv == 1 & ROSMAPIndiv$cogdx == 1, ]
ROSMAPIndivCtrls <- ROSMAPIndivCtrls[ROSMAPIndivCtrls$ceradsc == 4, ]
ROSMAPIndivCtrls <- ROSMAPIndivCtrls[ROSMAPIndivCtrls$braaksc <= 3, ]
ROSMAPIndivCtrls <- ROSMAPIndivCtrls[!is.na(ROSMAPIndivCtrls$cogdx), ]
dim(ROSMAPIndivCtrls)

# LBP
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn26337520

LBP <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/LBP"
LBP_files <- list.files(LBP, full.names = T)
LBPInd <- read.csv(LBP_files)

table(LBPInd$diagnosis)

# The only criteria here is diagnosis, so let's take controls and OCD and
# depression


isCtrl <- LBPInd$diagnosis == "control" | LBPInd$diagnosis == "obsessive compulsive disorder" | LBPInd$diagnosis == "Major Depressive Disorder"

LBPInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrl) + 1]

# Remove 89+
LBPInd <- LBPInd[LBPInd$ageDeath != "89+", ]
LBPInd$ageDeath <- as.numeric(LBPInd$ageDeath)
LBPInd$race <- gsub("Black or African American", "Black", LBPInd$race)

LBPIndCtrl <- LBPInd[LBPInd$diagn_4BrainClck == "Control", ]

write.csv(LBPInd, paste0(outDir, "LBP_ind_all.csv"))
write.csv(LBPIndCtrl, paste0(outDir, "LBP_ind_ctrls.csv"))

LBPCounts <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/expression/LBP/LBP_FlagshipPaper_featureCounts.csv")
LBPAssMeta <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/expression/LBP/LBP_assay_RNAseq_metadata.csv")
LBPBioMeta <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/expression/LBP/LBP_biospecimen_metadata.csv")

View(LBPIndCtrl)
View(LBPBioMeta)

LBPBioMeta$individualID[duplicated(LBPBioMeta$individualID)]

LBPIndCtrl$individualID %in% colnames(LBPCounts)

LBPCounts[1:10, 1:10]
# MIT ROSMAP Multiomics
################################################################################

# Same individuals as ROSMAP, so nothing new

# ROSMAP
################################################################################

# Already covered

# MSBB
################################################################################

# Already covered

# MayoRNAseq
################################################################################

# Already covered

# Jax.IU
################################################################################

Jax.IU <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/Jax.IU"
Jax.IU_files <- list.files(Jax.IU, full.names = T)

Jax.IUIndiv <- read.csv(Jax.IU_files)
table(Jax.IUIndiv$species)

# Is all mouse, so not useful

# ROSMAP Mamillary body
################################################################################

# Already covered

# FreshMicro
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study_Name=(FreshMicro)

freshMicro <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/freshMicro"
freshMicroFiles <- list.files(freshMicro, full.names = T)

freshMicroAss <- read.csv(freshMicroFiles[1])
freshMicroInd <- read.csv(freshMicroFiles[2])

# Use braak and ctrl diagnosis to filter
isCtrl <- freshMicroInd$diagnosis == "control" & freshMicroInd$Braak <= 3

freshMicroInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrl) + 1]

freshMicroInd <- freshMicroInd[!is.na(freshMicroInd$diagn_4BrainClck), ]

# Remove ambiguous ages and convert to numeric
freshMicroInd <- freshMicroInd[freshMicroInd$ageDeath != "90+", ]
freshMicroInd$ageDeath <- as.numeric(freshMicroInd$ageDeath)



# One control has as cause of deach complications of parkinson's. Remove it
freshMicroInd$diagn_4BrainClck[grepl("Parkinson",
                                     freshMicroInd$causeDeath)] <- "Rest"




freshMicroIndCtrl <- freshMicroInd[freshMicroInd$diagn_4BrainClck == "Control", ]




freshMicroCounts <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/expression/FreshMicro/RNAseq_count_matrix.csv")

# There is counts data of this dataset, so keep it

write.csv(freshMicroInd, file = paste0(outDir, "freshMicro_ind_all.csv"))
write.csv(freshMicroIndCtrl, file = paste0(outDir, "freshMicro_ind_ctrls.csv"))


# MayoPilotRNAseq
################################################################################

# Already covered

# SuperAgerEpiMap
################################################################################

# Built on MSBB cohort, so already covered

# MC-CAA
################################################################################

MC_CAA <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/MC-CAA"

MC_CAAFiles <- list.files(MC_CAA, full.names = T)


MC_CAAAss <- read.csv(MC_CAAFiles[1])
MC_CAAInd <- read.csv(MC_CAAFiles[2])

table(MC_CAAInd$diagnosis)

# Here is all AD. Might be used to asess brain age in disease individuals

# MayoHippocampus
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study_Name=(MayoHippocampus)

MayoHipp <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/MayoHippocampus"
MayoHippFiles <- list.files(MayoHipp, full.names = T)

MayoHippAss <- read.csv(MayoHippFiles[1])
MayoHippInd <- read.csv(MayoHippFiles[2])

MayoHippIndCtrl <- MayoHippInd[MayoHippInd$diagnosis == "control" & MayoHippInd$Braak <= 3, ]
MayoHippIndCtrl <- MayoHippIndCtrl[!is.na(MayoHippIndCtrl$diagnosis), ]
table(MayoHippIndCtrl$ageDeath)

# Data are FASTQs


MayoHarmIndCtrl$individualID
# RNAseq_Harmonization
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study_Name=(RNAseq_Harmonization)

# This is a harmonized version of ROSMAP, MSBB and Mayo, so samples should be 
# the same, but since they have been harmonized maybe is better to use these
# metadata.

RNAseqHarm <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/RNAseq_Harmonization"
RNAseqHarmFiles <- list.files(RNAseqHarm, full.names = T)
RNAseqHarmFiles <- RNAseqHarmFiles[grepl(".csv", RNAseqHarmFiles)]

MayoHarmInd <- read.csv(RNAseqHarmFiles[1])

# In the RNAseq Harmonization Mayo dataset there are several samples that
# don't have Braak score or Thal score, while in the original Mayo metadata
# they have. Replace them and filter based on that

MayoHarmInd$Braak[is.na(MayoHarmInd$Braak)] <- mayoIndiv$Braak[match(MayoHarmInd$individualID[is.na(MayoHarmInd$Braak)],
                                                                     mayoIndiv$individualID)]

MayoHarmInd$thal[is.na(MayoHarmInd$thal)] <- mayoIndiv$Thal[match(MayoHarmInd$individualID[is.na(MayoHarmInd$thal)],
                                                                     mayoIndiv$individualID)]


isCtrlMayo <- MayoHarmInd$diagnosis == "control" & MayoHarmInd$Braak <= 3 & MayoHarmInd$thal == 0

isCtrlMayo[is.na(isCtrlMayo)] <- F


MayoHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlMayo) + 1]

# Remove ambiguous ages and convert to numeric
MayoHarmInd <- MayoHarmInd[MayoHarmInd$ageDeath != "90+" & MayoHarmInd$ageDeath != "", ]

MayoHarmInd$ageDeath <- as.numeric(MayoHarmInd$ageDeath)

MayoHarmIndCtrl <- MayoHarmInd[MayoHarmInd$diagn_4BrainClck == "Control", ]

unique(MayoHarmIndCtrl$individualID)

# There are individualIDs duplicated, so check what they are
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][1], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][2], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][3], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][4], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][5], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][6], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][7], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][8], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][9], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][10], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][11], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][12], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][13], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][14], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][15], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][16], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][17], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][18], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][19], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][20], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][21], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][22], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][23], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][24], ]
MayoHarmIndCtrl[MayoHarmIndCtrl$individualID == MayoHarmIndCtrl$individualID[duplicated(MayoHarmIndCtrl$individualID)][25], ]

# The duplicated ones is because there is a cerebellum sample for that 
# individual. Probably we will need to remove it because it is a very different
# tissue, but we'll see.

MSBBHarmInd <- read.csv(RNAseqHarmFiles[2])

isCtrlMSBB <- MSBBHarmInd$Braak <= 3 & MSBBHarmInd$CERAD == 4
MSBBHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlMSBB) + 1]
MSBBHarmInd <- MSBBHarmInd[!is.na(MSBBHarmInd$diagn_4BrainClck), ]

# Remove ambiguous ages and convert to numeric
MSBBHarmInd <- MSBBHarmInd[MSBBHarmInd$ageDeath != "90+", ]
MSBBHarmInd$ageDeath <- as.numeric(MSBBHarmInd$ageDeath)

MSBBHarmIndCtrl <- MSBBHarmInd[MSBBHarmInd$diagn_4BrainClck == "Control", ]


ROSMAPHarmInd <- read.csv(RNAseqHarmFiles[3])

isCtrlROSMAP <- ROSMAPHarmInd$dcfdx_lv == 1 & ROSMAPHarmInd$cogdx == 1
isCtrlROSMAP <- isCtrlROSMAP & ROSMAPHarmInd$braaksc <= 3 & ROSMAPHarmInd$ceradsc == 4



ROSMAPHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlROSMAP) + 1]
ROSMAPHarmInd <- ROSMAPHarmInd[!is.na(ROSMAPHarmInd$diagn_4BrainClck), ]

# Remove ambiguous ages and, convert to numeric and round
ROSMAPHarmInd <- ROSMAPHarmInd[ROSMAPHarmInd$age_death != "90+", ]
ROSMAPHarmInd$age_death <- round(as.numeric(ROSMAPHarmInd$age_death))


ROSMAPHarmIndCtrl <- ROSMAPHarmInd[ROSMAPHarmInd$diagn_4BrainClck == "Control", ]


# Save the datasets in separated files.
write.csv(MayoHarmInd, file = paste0(outDir, "Mayo_Harm_parsed.csv"))
write.csv(MSBBHarmInd, file = paste0(outDir, "MSBB_Harm_parsed.csv"))
write.csv(ROSMAPHarmInd, file = paste0(outDir, "ROSMAP_Harm_parsed.csv"))


# Combine the three of them in a single dataframe

# Find the columns with common names and change the colnames with discrepancies
commCols <- colnames(MayoHarmIndCtrl)[colnames(MayoHarmInd) %in% colnames(MSBBHarmInd) & colnames(MayoHarmInd) %in% colnames(ROSMAPHarmInd)]

colnames(ROSMAPHarmInd)[!colnames(ROSMAPHarmInd) %in% commCols]

colnames(MayoHarmInd)[!colnames(MayoHarmInd) %in% commCols]

colnames(MSBBHarmInd)[!colnames(MayoHarmInd) %in% commCols]

colnames(ROSMAPHarmInd) <- gsub("apoe_genotype",
                                "apoeGenotype",
                                colnames(ROSMAPHarmInd))

ROSMAPHarmInd$msex <- c("female", "male")[ROSMAPHarmInd$msex + 1]
colnames(ROSMAPHarmInd) <- gsub("msex",
                                "sex",
                                colnames(ROSMAPHarmInd))

colnames(ROSMAPHarmInd) <- gsub("braaksc",
                                "Braak",
                                colnames(ROSMAPHarmInd))

colnames(ROSMAPHarmInd) <- gsub("age_death",
                                "ageDeath",
                                colnames(ROSMAPHarmInd))

ROSMAPHarmInd$race <- c("White",
                        "Black",
                        "Native_American",
                        "Native_Hawaiian",
                        "Asian",
                        "Other",
                        "Unknown")[ROSMAPHarmInd$race]

ROSMAPHarmInd$race[ROSMAPHarmInd$spanish == 1] <- "Hispanic"

MSBBHarmInd$race <- gsub("W", "White", MSBBHarmInd$race)
MSBBHarmInd$race <- gsub("B", "Black", MSBBHarmInd$race)
MSBBHarmInd$race <- gsub("H", "Hispanic", MSBBHarmInd$race)
MSBBHarmInd$race <- gsub("A", "Asian", MSBBHarmInd$race)
MSBBHarmInd$race <- gsub("U", "Unknown", MSBBHarmInd$race)

# Reassign commCols to account for the unifying changes and check that the
# format of the variables in the columns with common name is the same
commCols <- colnames(MayoHarmInd)[colnames(MayoHarmInd) %in% colnames(MSBBHarmInd) & colnames(MayoHarmInd) %in% colnames(ROSMAPHarmInd)]

head(ROSMAPHarmInd[, commCols])
head(MayoHarmInd[, commCols])
head(MSBBHarmInd[, commCols])

# Add a column of substudy and merge datasets.
ROSMAPHarmInd$substudy <- rep("ROSMAP", nrow(ROSMAPHarmInd))
MayoHarmInd$substudy <- rep("Mayo", nrow(MayoHarmInd))
MSBBHarmInd$substudy <- rep("MSBB", nrow(MSBBHarmInd))



commCols <- colnames(MayoHarmInd)[colnames(MayoHarmInd) %in%
                                          colnames(MSBBHarmInd) &
                                          colnames(MayoHarmInd) %in%
                                          colnames(ROSMAPHarmInd)]

RNAseqHarmDF <- rbind.data.frame(ROSMAPHarmInd[, commCols],
                                 MayoHarmInd[, commCols],
                                 MSBBHarmInd[, commCols])

write.csv(RNAseqHarmDF, file = paste0(outDir, "RNAseq_Harmonization_ind_all.csv"))

# Create a individual metadata file with both controls and patients



# RNAseqSampleSwap
################################################################################

# Built on samples from ROSMAP, MSBB and Mayo, so already covered

# SEA-AD
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study_Name=(SEA-AD)

SEA_AD <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/SEA-AD"
SEA_ADFiles <- list.files(SEA_AD, full.names = T)

SEA_ADInd <- read.csv(SEA_ADFiles[1])
SEA_ADAss <- read.csv(SEA_ADFiles[2])

isCtrl_SEA_AD <- SEA_ADInd$Consensus.clinical.diagnosis == "Control" & SEA_ADInd$Braak <= 3 & SEA_ADInd$CERAD == 0

SEA_ADInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrl_SEA_AD) + 1]

SEA_ADInd <- SEA_ADInd[SEA_ADInd$ageDeath != "90+", ]

SEA_ADInd$ageDeath <- as.numeric(SEA_ADInd$ageDeath)

SEA_ADIndCtrl <- SEA_ADInd[SEA_ADInd$diagn_4BrainClck == "Control", ]

table(SEA_ADIndCtrl$LATE.NC.stage)

SEA_ADIndCtrl

SEA_ADAss$donor %in% SEA_ADIndCtrl$individualID

# There are FPKM files of some samples of this dataset. Save ind metadata.

write.csv(SEA_ADInd, paste0(outDir, "SEA_AD_ind_all.csv"))
write.csv(SEA_ADIndCtrl, paste0(outDir, "SEA_AD_ind_ctrl.csv"))

# MCMPS
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study_Name=(MCMPS)

MCMPS <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/metadata/MCMPS"
MCMPSFiles <- list.files(MCMPS, full.names = T)
MCMPSAss <- read.csv(MCMPSFiles[1])
MCMPSInd <- read.csv(MCMPSFiles[2])

# Not useful, as there is no diagnosis and samples were obtained during brain
# surgeries

# BroadIPSCs
################################################################################

# Not useful, as they are iPSCs