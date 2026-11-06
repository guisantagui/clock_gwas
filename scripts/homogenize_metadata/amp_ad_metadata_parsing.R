################################################################################
# Brain clock: go through the individual metadata files of each one of the     #
# studies to determine which individuals are controls and which not and unify  #
# the information included in them.                                            #
################################################################################

setwd("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas")
outDir <- "./results/metadata_parsed/"

if (!dir.exists(outDir)){
    dir.create(outDir, recursive = T)
}

# LBP
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn26337520

LBP <- "./data/metadata/AMP_AD/LBP"
LBPInd <- read.csv(paste(LBP, "LBP_individual_metadata.csv", sep = "/"))

isCtrl <- LBPInd$diagnosis == "control"

LBPInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrl) + 1]

LBPInd$race <- gsub("Black or African American", "Black", LBPInd$race)

LBPSpc <- read.csv(file.path(LBP, "LBP_biospecimen_metadata.csv"))

LBP_WGS <- read.csv(file.path(LBP, "LBP_assay_WGS_metadata.csv"))

LBP_WGS$individualID <- LBPSpc$individualID[
        match(LBP_WGS$specimenID, LBPSpc$specimenID)
]

LBPInd$has_RNAseq <- LBPInd$individualID %in% LBPSpc$individualID[LBPSpc$assay == "rnaSeq"]
LBPInd$has_WGS <- LBPInd$individualID %in% LBP_WGS$individualID
LBPInd$has_snpArray <- LBPInd$individualID %in% LBPSpc$individualID[LBPSpc$assay == "snpArray"]

write.csv(LBPInd, paste0(outDir, "LBP_ind_all.csv"))

################################################################################
#                                                                              #
# RNAseq_Harmonization/WGS_Harmonization                                       #
#                                                                              #
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study_Name=(RNAseq_Harmonization)

# This is a harmonized version of ROSMAP, MSBB and Mayo, so samples should be 
# the same, but since they have been harmonized it's better to use these
# metadata.

RNAseqHarm <- "./data/metadata/AMP_AD/RNAseq_Harmonization"
WGSHarm <- "./data/metadata/AMP_AD/WGS_Harmonization/"

# Mayo Cohort
################################################################################

MayoHarmInd <- read.csv(paste(RNAseqHarm,
                              "/RNAseq_Harmonization_Mayo_combined_metadata.csv",
                              sep = "/"))

# In the RNAseq Harmonization Mayo dataset there are several samples that
# don't have Braak score or Thal score, while in the original Mayo metadata
# they have. Replace them and filter based on that.

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn9702085

#rnaSeqReprocessing <- "./data/metadata/rnaSeqReprocessing"
#rnaSeqReprocessing_files <- list.files(rnaSeqReprocessing, full.names = T)

mayoIndiv <- read.csv(paste(RNAseqHarm,
                            "MayoRNAseq_individual_metadata.csv",
                            sep = "/"))

mayoSpc <- read.csv(
        file.path(RNAseqHarm,
                  "MayoRNAseq_biospecimen_metadata.csv")
)

mayoRNAseqAssay <- read.csv(
        file.path(RNAseqHarm,
                  "MayoRNAseq_assay_RNAseq.csv")
)

sum(!MayoHarmInd$individualID %in% mayoIndiv$individualID)
sum(!mayoIndiv$individualID %in% MayoHarmInd$individualID)

# In mayo there are some individualIDs that are in the original Mayo metadata
# file that don't appear in the harmonized version. Assess if it's because
# these samples dont' have expression data

indIDs_with_genomics_and_rnaseq <- unique(
        intersect(
                mayoSpc$individualID[mayoSpc$assay == "snpArray" | mayoSpc$assay == "wholeGenomeSeq"], 
                mayoSpc$individualID[mayoSpc$assay == "rnaSeq"]
        )
)


sum(!indIDs_with_genomics_and_rnaseq %in% MayoHarmInd$individualID)
sum(!MayoHarmInd$individualID %in% indIDs_with_genomics_and_rnaseq)

# They do have expression data. Let's see if it's because of bad RIN

rnaseqharm_missing_inds_specID <- mayoSpc$specimenID[
        mayoSpc$individualID %in% indIDs_with_genomics_and_rnaseq[!indIDs_with_genomics_and_rnaseq %in% MayoHarmInd$individualID]
]

mayoRNAseqAssay[mayoRNAseqAssay$specimenID %in% rnaseqharm_missing_inds_specID, ]
# RIN is good

# So let's use mayoIndiv and add columns we are interested in
mayoIndiv_parsed <- mayoIndiv
mayoIndiv_parsed$ageDeath <- gsub("90_or_over", "90+", mayoIndiv_parsed$ageDeath)
mayoIndiv_parsed$diagn_4BrainClck <- mayoIndiv_parsed$diagnosis == "control" &
        mayoIndiv_parsed$Braak <= 3 & mayoIndiv_parsed$Thal == 0

mayoIndiv_parsed$diagn_4BrainClck[
        is.na(mayoIndiv_parsed$diagn_4BrainClck)
] <- FALSE

mayoIndiv_parsed$diagn_4BrainClck <- ifelse(
        mayoIndiv_parsed$diagn_4BrainClck, "Control", "Rest"
)


MayoHarmInd$Braak[is.na(MayoHarmInd$Braak)] <- mayoIndiv$Braak[match(MayoHarmInd$individualID[is.na(MayoHarmInd$Braak)],
                                                                     mayoIndiv$individualID)]

MayoHarmInd$thal[is.na(MayoHarmInd$thal)] <- mayoIndiv$Thal[match(MayoHarmInd$individualID[is.na(MayoHarmInd$thal)],
                                                                  mayoIndiv$individualID)]


isCtrlMayo <- MayoHarmInd$diagnosis == "control" & MayoHarmInd$Braak <= 3 & MayoHarmInd$thal == 0

isCtrlMayo[is.na(isCtrlMayo)] <- F



MayoHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlMayo) + 1]

# RNAseq harmonization and WGS harmonization have been run through the same
# pipelines, so let's focus on what's in there

# Load WGS assay data 
MayoHarm_wgs_assay <- read.csv(file.path(WGSHarm, "MayoRNAseq_assay_wholeGenomeSeq_metadata.csv"))


mayoIndiv_parsed <- mayoIndiv_parsed[
        mayoIndiv_parsed$individualID %in% MayoHarmInd$individualID |
        mayoIndiv_parsed$individualID %in% mayoSpc$individualID[match(MayoHarm_wgs_assay$specimenID, mayoSpc$specimenID)],
]

mayoIndiv_parsed$has_RNAseq <- mayoIndiv_parsed$individualID %in% MayoHarmInd$individualID
mayoIndiv_parsed$has_WGS <- mayoIndiv_parsed$individualID %in% mayoSpc$individualID[match(MayoHarm_wgs_assay$specimenID, mayoSpc$specimenID)]
mayoIndiv_parsed$has_snpArray <- mayoIndiv_parsed$individualID %in% mayoSpc$individualID[mayoSpc$assay == "snpArray"]

#mayoIndiv_parsed$has_RNAseq <- mayoIndiv_parsed$individualID %in% mayoSpc$individualID[mayoSpc$assay == "rnaSeq"]
#mayoIndiv_parsed$has_WGS <- mayoIndiv_parsed$individualID %in% mayoSpc$individualID[mayoSpc$assay == "wholeGenomeSeq"]
#mayoIndiv_parsed$has_snpArray <- mayoIndiv_parsed$individualID %in% mayoSpc$individualID[mayoSpc$assay == "snpArray"]

# MSBB Cohort
################################################################################

MSBBHarmInd <- read.csv(paste(RNAseqHarm,
                              "/RNAseq_Harmonization_MSBB_combined_metadata.csv",
                              sep = "/"))

isCtrlMSBB <- MSBBHarmInd$Braak <= 3 & MSBBHarmInd$CERAD == 4
MSBBHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlMSBB) + 1]
MSBBHarmInd$diagn_4BrainClck[is.na(MSBBHarmInd$diagn_4BrainClck)] <- "Rest"

msbb_dir <- "./data/metadata/AMP_AD/MSBB"
msbbIndiv <- read.csv(file.path(msbb_dir, "MSBB_individual_metadata_harmonized.csv"))

msbbSpc <- read.csv(
        file.path(msbb_dir,
                  "MSBB_biospecimen_metadata.csv")
)

msbbRNAseqAssay <- read.csv(
        file.path(msbb_dir,
                  "MSBB_assay_rnaSeq_metadata.csv")
)

sum(!unique(MSBBHarmInd$individualID) %in% msbbIndiv$individualID)
sum(!msbbIndiv$individualID %in% MSBBHarmInd$individualID)

# The harmonized MSBB has 1 individuals that's not in MSBB's original
# individual file
# file that don't appear in the harmonized version. Assess if it's because
# The original MSBB individual file has 62 individuals not in harmonized MSBB

indIDs_with_genomics_and_rnaseq <- unique(
        intersect(
                msbbSpc$individualID[msbbSpc$assay == "snpArray" | msbbSpc$assay == "wholeGenomeSeq"], 
                msbbSpc$individualID[msbbSpc$assay == "rnaSeq"]
        )
)


sum(!indIDs_with_genomics_and_rnaseq %in% MSBBHarmInd$individualID)
sum(!unique(MSBBHarmInd$individualID) %in% indIDs_with_genomics_and_rnaseq)

# All the individuals with both RNAseq and genomics are ion the harmonized
# individual file. 19 of the harmonized individuals that are not in the
# individuals from the original fiule with both genomics and transcriptomics

# There are 69 samples 

# So let's use MSBB as base, and add columns of interest and verify if these
# individuals have genomics data as per the WGS harmonization study
msbbIndiv_parsed <- msbbIndiv
msbbIndiv_parsed$diagn_4BrainClck <- msbbIndiv_parsed$Braak <= 3 & msbbIndiv_parsed$CERAD == 4

mayoIndiv_parsed$diagn_4BrainClck[
        is.na(mayoIndiv_parsed$diagn_4BrainClck)
] <- FALSE

msbbHarm_ctrls <- unique(MSBBHarmInd$individualID[MSBBHarmInd$diagn_4BrainClck == "Control"])

msbbHarm_discrepant_controls <- msbbHarm_ctrls[msbbHarm_ctrls %in% msbbIndiv_parsed$individualID[!msbbIndiv_parsed$diagn_4BrainClck]]


MSBBHarmInd[MSBBHarmInd$individualID %in% msbbHarm_discrepant_controls, ] %>% distinct(individualID, CERAD, Braak)

msbbIndiv_parsed[msbbIndiv_parsed$individualID %in% msbbHarm_discrepant_controls, ]

msbb_harm_clinical <- MSBBHarmInd %>% distinct(individualID, CERAD, Braak, ageDeath)
msbb_orig_clinical <- msbbIndiv[match(msbb_harm_clincial$individualID, msbbIndiv$individualID), c("individualID", "CERAD", "Braak", "ageDeath")]

colnames(msbb_harm_clinical) <- paste0(colnames(msbb_harm_clinical), "_harm")
colnames(msbb_orig_clinical) <- paste0(colnames(msbb_orig_clinical), "_orig")

msbb_harm_clinical$ageDeath_harm <- gsub("\\..*", "", msbb_harm_clinical$ageDeath_harm)

msbb_clinical_comb <- cbind.data.frame(msbb_harm_clinical, msbb_orig_clinical)
msbb_clinical_comb[(msbb_clinical_comb$Braak_harm != msbb_clinical_comb$Braak_orig) | (msbb_clinical_comb$CERAD_harm != msbb_clinical_comb$CERAD_orig) | (msbb_clinical_comb$ageDeath_harm != msbb_clinical_comb$ageDeath_orig), ]
msbbIndiv[match(msbb_harm_clincial$individualID, msbbIndiv$individualID), ]$Braak == msbb_harm_clinical$Braak


mayoIndiv_parsed$diagn_4BrainClck <- ifelse(
        mayoIndiv_parsed$diagn_4BrainClck, "Control", "Rest"
)


MayoHarmInd$Braak[is.na(MayoHarmInd$Braak)] <- mayoIndiv$Braak[match(MayoHarmInd$individualID[is.na(MayoHarmInd$Braak)],
                                                                     mayoIndiv$individualID)]

MayoHarmInd$thal[is.na(MayoHarmInd$thal)] <- mayoIndiv$Thal[match(MayoHarmInd$individualID[is.na(MayoHarmInd$thal)],
                                                                  mayoIndiv$individualID)]


isCtrlMayo <- MayoHarmInd$diagnosis == "control" & MayoHarmInd$Braak <= 3 & MayoHarmInd$thal == 0

isCtrlMayo[is.na(isCtrlMayo)] <- F



MayoHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlMayo) + 1]

# RNAseq harmonization and WGS harmonization have been run through the same
# pipelines, so let's focus on what's in there

# Load WGS assay data 
MayoHarm_wgs_assay <- read.csv(file.path(WGSHarm, "MayoRNAseq_assay_wholeGenomeSeq_metadata.csv"))


mayoIndiv_parsed <- mayoIndiv_parsed[
        mayoIndiv_parsed$individualID %in% MayoHarmInd$individualID |
                mayoIndiv_parsed$individualID %in% mayoSpc$individualID[match(MayoHarm_wgs_assay$specimenID, mayoSpc$specimenID)],
]

mayoIndiv_parsed$has_RNAseq <- mayoIndiv_parsed$individualID %in% MayoHarmInd$individualID
mayoIndiv_parsed$has_WGS <- mayoIndiv_parsed$individualID %in% mayoSpc$individualID[match(MayoHarm_wgs_assay$specimenID, mayoSpc$specimenID)]
mayoIndiv_parsed$has_snpArray <- mayoIndiv_parsed$individualID %in% mayoSpc$individualID[mayoSpc$assay == "snpArray"]














# Remove ambiguous ages and convert to numeric
#MSBBHarmInd <- MSBBHarmInd[MSBBHarmInd$ageDeath != "90+", ]
#MSBBHarmInd$ageDeath <- as.numeric(MSBBHarmInd$ageDeath)

#MSBBHarmIndCtrl <- MSBBHarmInd[MSBBHarmInd$diagn_4BrainClck == "Control", ]


ROSMAPHarmInd <- read.csv(paste(RNAseqHarm,
                              "/RNAseq_Harmonization_ROSMAP_combined_metadata.csv",
                              sep = "/"))

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