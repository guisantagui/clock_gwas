################################################################################
# GWAS of Brain clock: go through the individual metadata files of each one of #
# the studies to determine which individuals are controls and which not, and   #
# generate a unified individual metadata file that we can use to assess what   # 
# individuals have both RNAseq and genomics data                               #
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
sum(!unique(MayoHarmInd$individualID) %in% indIDs_with_genomics_and_rnaseq)


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

MayoHarmInd$cohort <- "MayoClinic"
mayoIndiv_parsed <- MayoHarmInd %>% distinct(individualID, cohort, sex, race, ageDeath, apoeGenotype, pmi, diagnosis, Braak, thal, diagn_4BrainClck)
mayoIndiv_parsed$has_RNAseq <- mayoIndiv_parsed$individualID %in% MayoHarmInd$individualID
mayoIndiv_parsed$has_WGS <- mayoIndiv_parsed$individualID %in% mayoSpc$individualID[match(MayoHarm_wgs_assay$specimenID, mayoSpc$specimenID)]
mayoIndiv_parsed$has_snpArray <- mayoIndiv_parsed$individualID %in% mayoSpc$individualID[mayoSpc$assay == "snpArray"]

# MSBB Cohort
################################################################################
msbb_dir <- "./data/metadata/AMP_AD/MSBB"

MSBBHarmInd <- read.csv(paste(RNAseqHarm,
                              "/RNAseq_Harmonization_MSBB_combined_metadata.csv",
                              sep = "/"))

# In the RNAseq Harmonization Mayo dataset there are several samples that
# don't have Braak score or Thal score, while in the original Mayo metadata
# they have. Replace them and filter based on that.

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn9702085

#rnaSeqReprocessing <- "./data/metadata/rnaSeqReprocessing"
#rnaSeqReprocessing_files <- list.files(rnaSeqReprocessing, full.names = T)

msbbIndiv <- read.csv(file.path(msbb_dir, "MSBB_individual_metadata.csv"))

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

# In mayo there are some individualIDs that are in the original Mayo metadata
# file that don't appear in the harmonized version. Assess if it's because
# these samples dont' have expression data

indIDs_with_genomics_and_rnaseq <- unique(
        intersect(
                msbbSpc$individualID[msbbSpc$assay == "snpArray" | msbbSpc$assay == "wholeGenomeSeq"], 
                msbbSpc$individualID[msbbSpc$assay == "rnaSeq"]
        )
)


sum(!indIDs_with_genomics_and_rnaseq %in% MSBBHarmInd$individualID)
sum(!unique(MSBBHarmInd$individualID) %in% indIDs_with_genomics_and_rnaseq)

# In the brain clock we fucked up and used the same CERAD scale as in rosmap
# So individuals used in training were according to CERAD "Possible AD". Still,
# that would fall in the No AD classification, so as per the Braak index they
# still categorize as controls. But let's fix CERAD so its the same in ROSMAP
# and MSBB

isCtrlMSBB <- MSBBHarmInd$Braak <= 3 & MSBBHarmInd$CERAD == 4
MSBBHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlMSBB) + 1]
#MSBBHarmInd <- MSBBHarmInd[!is.na(MSBBHarmInd$diagn_4BrainClck), ]
MSBBHarmInd$diagn_4BrainClck[is.na(MSBBHarmInd$diagn_4BrainClck)] <- "Rest"

table(MSBBHarmInd$CERAD[MSBBHarmInd$Braak <= 1])

CERAD_dict <- data.frame(
        cerad = c("Definite AD", "Probable AD", "Possible AD", "No AD"),
        rosmap_code = 1:4,
        msbb_code = c(2, 3, 4, 1)
)

# Also, apparently there was a problem detected in August 2025 with 
# MSBB metadata, so let's use the specific one from the study, which has been
# updated. It affected CERAD, Braak, sex, APOE genotype, CDR and PMI

MSBBHarmInd$ageDeath <- msbbIndiv$ageDeath[
        match(MSBBHarmInd$individualID,
              msbbIndiv$individualID)
]

MSBBHarmInd$CERAD <- msbbIndiv$CERAD[
        match(MSBBHarmInd$individualID,
              msbbIndiv$individualID)
]

MSBBHarmInd$Braak <- msbbIndiv$Braak[
        match(MSBBHarmInd$individualID,
              msbbIndiv$individualID)
]

MSBBHarmInd$sex <- msbbIndiv$sex[
        match(MSBBHarmInd$individualID,
              msbbIndiv$individualID)
]

MSBBHarmInd$apoeGenotype <- msbbIndiv$apoeGenotype[
        match(MSBBHarmInd$individualID,
              msbbIndiv$individualID)
]

MSBBHarmInd$CDR <- msbbIndiv$CDR[
        match(MSBBHarmInd$individualID,
              msbbIndiv$individualID)
]

MSBBHarmInd$pmi <- msbbIndiv$pmi[
        match(MSBBHarmInd$individualID,
              msbbIndiv$individualID)
]

MSBBHarmInd$CERAD <- CERAD_dict$cerad[
        match(MSBBHarmInd$CERAD, CERAD_dict$msbb_code)
]

# RNAseq harmonization and WGS harmonization have been run through the same
# pipelines, so let's focus on what's in there

# Load WGS assay data 
MSBBHarm_wgs_assay <- read.csv(file.path(WGSHarm, "MSBB_assay_wholeGenomeSeq_metadata.csv"))

MSBBHarmInd$cohort <- "MSBB"
msbbIndiv_parsed <- MSBBHarmInd %>% distinct(individualID, cohort, sex, race, ageDeath, apoeGenotype, pmi, Braak, CERAD, CDR, plaqueMean, diagn_4BrainClck)
msbbIndiv_parsed <- msbbIndiv_parsed[msbbIndiv_parsed$individualID != "", ]
msbbIndiv_parsed$has_RNAseq <- msbbIndiv_parsed$individualID %in% MSBBHarmInd$individualID
msbbIndiv_parsed$has_WGS <- msbbIndiv_parsed$individualID %in% msbbSpc$individualID[match(MSBBHarm_wgs_assay$specimenID, msbbSpc$specimenID)]
msbbIndiv_parsed$has_snpArray <- msbbIndiv_parsed$individualID %in% msbbSpc$individualID[msbbSpc$assay == "snpArray"]

# ROSMAP Cohort
################################################################################
rosmap_dir <- "./data/metadata/AMP_AD/ROSMAP/"

ROSMAPHarmInd <- read.csv(paste(RNAseqHarm,
                                "/RNAseq_Harmonization_ROSMAP_combined_metadata.csv",
                                sep = "/"))

isCtrlROSMAP <- ROSMAPHarmInd$dcfdx_lv == 1 & ROSMAPHarmInd$cogdx == 1
isCtrlROSMAP <- isCtrlROSMAP & ROSMAPHarmInd$braaksc <= 3 & ROSMAPHarmInd$ceradsc == 4



ROSMAPHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlROSMAP) + 1]

ROSMAPHarmInd$age_death[ROSMAPHarmInd$age_death == ""] <- NA

ROSMAPHarmInd$age_death[
        ROSMAPHarmInd$age_death != "90+" & !is.na(ROSMAPHarmInd$age_death)
] <- as.character(
        round(
                as.numeric(
                        ROSMAPHarmInd$age_death[
                                ROSMAPHarmInd$age_death != "90+" &
                                !is.na(ROSMAPHarmInd$age_death)
                        ]
                )
        )
)

ROSMAPHarm_wgs_assay <- read.csv(file.path(WGSHarm, "ROSMAP_assay_wholeGenomeSeq_metadata.csv"))
rosmapSpc <- read.csv(file.path(rosmap_dir, "ROSMAP_biospecimen_metadata.csv"))

rosmapIndiv_parsed <- ROSMAPHarmInd %>% distinct(individualID, Study, msex, race, spanish, age_death, apoe_genotype, pmi, braaksc, ceradsc, cogdx, dcfdx_lv, diagn_4BrainClck)
rosmapIndiv_parsed <- rosmapIndiv_parsed[rosmapIndiv_parsed$individualID != "", ]
rosmapIndiv_parsed$has_RNAseq <- rosmapIndiv_parsed$individualID %in% ROSMAPHarmInd$individualID
rosmapIndiv_parsed$has_WGS <- rosmapIndiv_parsed$individualID %in% rosmapSpc$individualID[match(ROSMAPHarm_wgs_assay$specimenID, rosmapSpc$specimenID)]
rosmapIndiv_parsed$has_snpArray <- rosmapIndiv_parsed$individualID %in% rosmapSpc$individualID[rosmapSpc$assay == "snpArray"]

rosmapIndiv_parsed$ceradsc <- CERAD_dict$cerad[
        match(rosmapIndiv_parsed$ceradsc, CERAD_dict$rosmap_code)
]

################################################################################
#                                                                              #
# Merge all the metadata files                                                 #
#                                                                              #
################################################################################

commCols <- colnames(rosmapIndiv_parsed) %>%
        intersect(colnames(mayoIndiv_parsed)) %>%
        intersect(colnames(msbbIndiv_parsed))

colnames(mayoIndiv_parsed)[!colnames(mayoIndiv_parsed) %in% commCols]
head(mayoIndiv_parsed[, colnames(mayoIndiv_parsed)[!colnames(mayoIndiv_parsed) %in% commCols]])

colnames(rosmapIndiv_parsed) <- gsub(
        "apoe_genotype",
        "apoeGenotype",
        colnames(rosmapIndiv_parsed)
)

rosmapIndiv_parsed$msex <- ifelse(rosmapIndiv_parsed$msex == 0, "female", "male")

colnames(rosmapIndiv_parsed) <- gsub(
        "msex",
        "sex",
        colnames(rosmapIndiv_parsed)
)

colnames(rosmapIndiv_parsed) <- gsub(
        "ceradsc",
        "CERAD",
        colnames(rosmapIndiv_parsed)
)

colnames(rosmapIndiv_parsed) <- gsub(
        "braaksc",
        "Braak",
        colnames(rosmapIndiv_parsed)
)

colnames(rosmapIndiv_parsed) <- gsub(
        "age_death",
        "ageDeath",
        colnames(rosmapIndiv_parsed)
)

rosmapIndiv_parsed$race <- c(
        "White",
        "Black",
        "Native_American",
        "Native_Hawaiian",
        "Asian",
        "Other",
        "Unknown")[rosmapIndiv_parsed$race]

rosmapIndiv_parsed$race[rosmapIndiv_parsed$spanish == 1] <- "Hispanic"

dcfdx_lv_meaning <- data.frame(value = 1:6,
                               diagnosis = c(
                                       "control",
                                       "mild cognitive impairment",
                                       "mild cognitive impairment",
                                       "Alzheimer Disease",
                                       "Alzheimer Disease",
                                       "other dementia"))
rosmapIndiv_parsed$diagnosis <- dcfdx_lv_meaning$diagnosis[
        match(rosmapIndiv_parsed$dcfdx_lv, dcfdx_lv_meaning$value)
]

msbb_race_code <- data.frame(
        code = c("A", "B", "H", "O", "U", "W"),
        race = c("Asian", "Black", "Hispanic", "Other", "Unknown", "White")
)
msbbIndiv_parsed$race <- msbb_race_code$race[match(msbbIndiv_parsed$race, msbb_race_code$code)]

# in MSBB PMI is in minutes. Pass to hours
msbbIndiv_parsed$pmi <- msbbIndiv_parsed$pmi / 60

colnames(rosmapIndiv_parsed) <- gsub("Study", "cohort", colnames(rosmapIndiv_parsed))

head(msbbIndiv_parsed)
head(rosmapIndiv_parsed)

rosmapIndiv_parsed$CDR <- NA
rosmapIndiv_parsed$plaqueMean <- NA
rosmapIndiv_parsed$thal <- NA

rosmapIndiv_parsed <- rosmapIndiv_parsed[,
        c(
                "individualID",
                "cohort",
                "sex",
                "race",
                "ageDeath",
                "apoeGenotype",
                "pmi",
                "diagnosis",
                "Braak",
                "CERAD",
                "CDR",
                "plaqueMean",
                "thal",
                "diagn_4BrainClck",
                "has_RNAseq",
                "has_WGS",
                "has_snpArray"
        )
]

msbbIndiv_parsed$thal <- NA

msbbIndiv_parsed$diagnosis <- ifelse(
        msbbIndiv_parsed$CERAD == "No AD" & msbbIndiv_parsed$Braak <= 3,
        "control",
        "Alzheimer Disease"
)

head(msbbIndiv_parsed)
head(rosmapIndiv_parsed)

msbbIndiv_parsed <- msbbIndiv_parsed[, colnames(rosmapIndiv_parsed)]

mayoIndiv_parsed$CERAD <- NA
mayoIndiv_parsed$CDR <- NA
mayoIndiv_parsed$plaqueMean <- NA
mayoIndiv_parsed <- mayoIndiv_parsed[, colnames(rosmapIndiv_parsed)]

head(LBPInd)

table(LBPInd$PMI)

lbp_pmi_split <- sapply(
        LBPInd$PMI, function(x) strsplit(x, split = ":")
)

lbp_pmi_length <- sapply(
        LBPInd$PMI, function(x) length(strsplit(x, split = ":")[[1]])
)

parsed_pmis <- numeric(length = length(LBPInd$PMI))
for(i in seq_along(LBPInd$PMI)){
        pmi <- LBPInd$PMI[i]
        if (pmi == ""){
                pmi <- NA
        }else{
                pmi_split <- strsplit(pmi, split = ":")[[1]]
                if (length(pmi_split) > 1){
                        pmi <- 0
                        for (j in seq_along(pmi_split)){
                                pmi <- pmi + as.numeric(pmi_split[j]) * (1 / (60^(j - 1)))
                        }
                }else{
                        pmi <- as.numeric(pmi_split[1]) * 24
                }
                
        }
        parsed_pmis[i] <- pmi
}

head(LBPInd)

LBPInd$PMI <- parsed_pmis
table(LBPInd$diagnosis)

head(LBPInd)

head(mayoIndiv_parsed)
head(LBPInd)

colnames(LBPInd) <- gsub(
        "PMI",
        "pmi",
        colnames(LBPInd)
)
LBPInd$CDR <- NA
LBPInd$plaqueMean <- NA
LBPInd$thal <- NA

LBP_Ind_parsed <- LBPInd[, colnames(rosmapIndiv_parsed)]

LBP_Ind_parsed$substudy <- "LBP"
rosmapIndiv_parsed$substudy <- "ROSMAP"
mayoIndiv_parsed$substudy <- "Mayo"
msbbIndiv_parsed$substudy <- "MSBB"

amp_ad_ind_mdat_merged <- rbind.data.frame(
        LBP_Ind_parsed,
        rosmapIndiv_parsed,
        mayoIndiv_parsed,
        msbbIndiv_parsed
)

sum(amp_ad_ind_mdat_merged$has_RNAseq & (amp_ad_ind_mdat_merged$has_WGS | amp_ad_ind_mdat_merged$has_snpArray))

clock_train_test_f <- "./data/metadata/brain_clock/train_test.csv"

clock_train_test <- read.csv(clock_train_test_f, row.names = 1)

clock_training_inds <- unique(
        clock_train_test$individualID[clock_train_test$train_test_split == "train"]
)
table(clock_train_test$substudy[match(clock_training_inds, clock_train_test$individualID)])

amp_ad_ind_mdat_merged$used_in_clock_training <- make.names(
        amp_ad_ind_mdat_merged$individualID
) %in% make.names(clock_training_inds)


amp_ad_ind_mdat_merged <- amp_ad_ind_mdat_merged[, c(1:2, 18, 3:17, 19)]

write.csv(
        amp_ad_ind_mdat_merged,
        file = file.path(outDir, "amp_ad_individual_mdat.csv")
)

sum(amp_ad_ind_mdat_merged$has_RNAseq & (amp_ad_ind_mdat_merged$has_WGS | amp_ad_ind_mdat_merged$has_snpArray) & !amp_ad_ind_mdat_merged$used_in_clock_training)



head(rosmapSpc[rosmapSpc$assay == "wholeGenomeSeq", ]$tissue)
rosmapSpc[rosmapSpc$assay == "wholeGenomeSeq", ]$individualID[duplicated(rosmapSpc[rosmapSpc$assay == "wholeGenomeSeq", ]$individualID)]


rosmap_wgs_biosp <- rosmapSpc[rosmapSpc$assay == "wholeGenomeSeq", ]
head(rosmap_wgs_biosp)

head(ROSMAPHarm_wgs_assay)
head(rosmapSpc)


head(MSBBHarm_wgs_assay)
head(msbbSpc)

head(MayoHarm_wgs_assay)
head(mayoSpc)


MayoHarm_wgs_assay$individualID <- mayoSpc$individualID[match(MayoHarm_wgs_assay$specimenID, mayoSpc$specimenID)]
MSBBHarm_wgs_assay$individualID <- msbbSpc$individualID[match(MSBBHarm_wgs_assay$specimenID, msbbSpc$specimenID)]
ROSMAPHarm_wgs_assay$individualID <- rosmapSpc$individualID[match(ROSMAPHarm_wgs_assay$specimenID, rosmapSpc$specimenID)]

head(MayoHarm_wgs_assay)
head(MSBBHarm_wgs_assay)
head(ROSMAPHarm_wgs_assay)

ROSMAPHarm_wgs_assay$Id <- NA
MSBBHarm_wgs_assay$Id <- NA
MayoHarm_wgs_assay$Id <- NA

head(LBP_WGS)

commCols_wgsMdat <- colnames(MayoHarm_wgs_assay) %>%
        intersect(colnames(MSBBHarm_wgs_assay)) %>%
        intersect(colnames(ROSMAPHarm_wgs_assay)) %>%
        intersect(colnames(LBP_WGS))

commCols_wgsMdat <- commCols_wgsMdat[c(1:2, 13:14, 3:12)]

wgs_amp_ad_mdat_comb <- rbind.data.frame(
        MayoHarm_wgs_assay[, commCols_wgsMdat],
        MSBBHarm_wgs_assay[, commCols_wgsMdat],
        ROSMAPHarm_wgs_assay[, commCols_wgsMdat],
        LBP_WGS[, commCols_wgsMdat]
)

write.csv(
        wgs_amp_ad_mdat_comb,
        file = file.path(outDir, "wgs_amp_ad_mdat_comb.csv")
)
